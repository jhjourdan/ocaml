#include <math.h>
#include <string.h>
#include "caml/memprof.h"
#include "caml/backtrace.h"
#include "caml/signals.h"
#include "caml/stacks.h"
#include "caml/fail.h"
#include "caml/memory.h"
#include "caml/alloc.h"
#include "caml/hash.h"
#include "caml/callback.h"
#ifdef NATIVE_CODE
#include "stack.h"
#endif

static uint32_t mt_state[624];
static uint32_t mt_index;

/* [lambda] is the mean number of samples for each allocated word (including
   block headers. */
static double lambda = 0;
static double lambda_rec = INFINITY;

static double next_sample_young;
value* caml_memprof_young_limit;

struct tracked_block {
  value block;
#ifdef NATIVE_CODE
  uint32_t alloc_frame_pos;
#endif
  uint32_t occurences;
  value callstack;
};

struct tracked_block* tracked_blocks = NULL;
uintnat tracked_blocks_end = 0;
static uintnat size = 0, old = 0;

/* Taken from :
   https://github.com/jhjourdan/SIMD-math-prims
*/
inline float logapprox(float val) {
  union { float f; int i; } valu;
  float exp, addcst, x;
  valu.f = val;
  exp = valu.i >> 23;
  addcst = val > 0 ? -89.970756366f : -(float)INFINITY;
  valu.i = (valu.i & 0x7FFFFF) | 0x3F800000;
  x = valu.f;

  return
    x * (3.529304993f + x * (-2.461222105f +
      x * (1.130626167f + x * (-0.288739945f +
        x * 3.110401639e-2f))))
    + (addcst + 0.69314718055995f*exp);
}

inline float expapprox(float val) {
  union { int i; float f; } xu, xu2;
  float val2, val3, val4, b;
  int val4i;
  val2 = 12102203.1615614f*val+1065353216.f;
  val3 = val2 < 2139095040.f ? val2 : 2139095040.f;
  val4 = val3 > 0.f ? val3 : 0.f;
  val4i = (int) val4;
  xu.i = val4i & 0x7F800000;
  xu2.i = (val4i & 0x7FFFFF) | 0x3F800000;
  b = xu2.f;

  return
    xu.f * (0.510397365625862338668154f + b *
            (0.310670891004095530771135f + b *
             (0.168143436463395944830000f + b *
              (-2.88093587581985443087955e-3f + b *
               1.3671023382430374383648148e-2f))));
}

static double mt_generate_uniform(void) {
  int i;
  uint32_t y;

  /* Mersenne twister PRNG */
  if (mt_index == 624) {
    for(i = 0; i < 227; i++) {
      y = (mt_state[i] & 0x80000000) + (mt_state[i+1] & 0x7fffffff);
      mt_state[i] = mt_state[i+397] ^ (y >> 1) ^ ((-(y&1)) & 0x9908b0df);
    }
    for(i = 227; i < 623; i++) {
      y = (mt_state[i] & 0x80000000) + (mt_state[i+1] & 0x7fffffff);
      mt_state[i] = mt_state[i-227] ^ (y >> 1) ^ ((-(y&1)) & 0x9908b0df);
    }
    y = (mt_state[623] & 0x80000000) + (mt_state[0] & 0x7fffffff);
    mt_state[623] = mt_state[396] ^ (y >> 1) ^ ((-(y&1)) & 0x9908b0df);
    mt_index = 0;
  }

  y = mt_state[mt_index];
  y = y ^ (y >> 11);
  y = y ^ ((y << 7) & 0x9d2c5680);
  y = y ^ ((y << 15) & 0xefc60000);
  y = y ^ (y >> 18);

  mt_index++;
  return y*2.3283064365386962890625e-10 + /* 2^-32 */
          1.16415321826934814453125e-10; /* 2^-33 */
}

static double mt_generate_exponential() {
  Assert(lambda >= 0 && !isinf(lambda));

  if(lambda == 0)
    return INFINITY;

  double res = -logapprox(mt_generate_uniform()) * lambda_rec;
  if(res < 0) return 0;
  return res;
}

/* Max returned value : 2^30-2 */
static uint32_t mt_generate_poisson(double lambda) {
  Assert(lambda >= 0 && !isinf(lambda));

  if(lambda == 0)
    return 0;

  if(lambda < 20) {
    double p;
    uint32_t k;
    k = 0;
    p = expapprox(lambda);
    do {
      k++;
      p *= mt_generate_uniform();
    } while(p > 1);
    return k-1;
  } else {
    double c, beta, alpha, k;
    c = 0.767 - 3.36/lambda;
    beta = 1./sqrt((3./(M_PI*M_PI))*lambda);
    alpha = beta*lambda;
    k = logapprox(c) - lambda - logapprox(beta);
    while(1) {
      double u, x, n, v, y, y2;
      u = mt_generate_uniform();
      x = (alpha - logapprox((1.0 - u)/u))/beta;
      n = floor(x + 0.5);
      if(n < 0.)
        continue;

      v = mt_generate_uniform();
      y = alpha - beta*x;
      y2 = 1. + expapprox(y);
      if(y + logapprox(v/(y2*y2)) < k + n*logapprox(lambda) - lgammaf(n+1))
        return n > ((1<<30)-2) ? ((1<<30)-2) : n;
    }
  }
}

static void update_limit(void) {
  Assert(next_sample_young >= caml_young_alloc_end - caml_young_ptr);
  if(next_sample_young > caml_young_alloc_end - caml_young_alloc_start + Max_young_whsize)
    caml_memprof_young_limit = caml_young_alloc_start - Max_young_wosize;
  else
    caml_memprof_young_limit = caml_young_alloc_end - (uintnat)next_sample_young;

  caml_young_limit =
    caml_memprof_young_limit < caml_young_limit ?
    caml_young_limit : caml_memprof_young_limit;
}

static void renew_sample(void) {
  next_sample_young =
    caml_young_alloc_end - caml_young_ptr + mt_generate_exponential();

  update_limit();
}

void caml_memprof_reinit(void) {
  int i;

  if(tracked_blocks == NULL) {
    mt_index = 624;
    mt_state[0] = 42;
    for(i = 1; i < 624; i++)
      mt_state[i] = 0x6c078965 * (mt_state[i-1] ^ (mt_state[i-1] >> 30)) + i;

    tracked_blocks =
      (struct tracked_block*)malloc(1024 * sizeof(struct tracked_block));
    if(tracked_blocks == NULL)
      caml_fatal_error("out of memory");
    size = 1024;
  }

  renew_sample();
}

static void push(const struct tracked_block* blk) {
  if(tracked_blocks_end >= size) {
    struct tracked_block* p =
      (struct tracked_block*)realloc(tracked_blocks,
        size * 2 * sizeof(struct tracked_block));

    if(p == NULL)
      return;

    tracked_blocks = p;

    size *= 2;
  }

  Assert(tracked_blocks_end < size);
  Assert(Is_in_heap_or_young(blk->block));
  tracked_blocks[tracked_blocks_end++] = *blk;
}

static void shrink(void) {
  uintnat newsize = size;
  while(newsize >= 2048 && tracked_blocks_end * 8 < newsize)
    newsize /= 2;
  if(newsize != size) {
    struct tracked_block* p = (struct tracked_block*)realloc(tracked_blocks,
                                        newsize * sizeof(struct tracked_block));
    if(p != NULL) {
      tracked_blocks = p;
      size = newsize;
    }
  }
}

static value memprof_callback;

CAMLprim value caml_memprof_register_callback(value callback) {
  CAMLparam1(callback);
  memprof_callback = callback;
  caml_register_global_root(&memprof_callback);
  CAMLreturn(Val_unit);
}

void caml_memprof_track_young(uintnat wosize) {
  uintnat whsize = Whsize_wosize(wosize);
  CAMLparam0();
  CAMLlocal1(callstack);
  double rest = caml_young_alloc_end - caml_young_ptr - next_sample_young;
  struct tracked_block blk;

  Assert(rest > 0);

  caml_young_ptr += whsize;

  renew_sample();

#ifdef NATIVE_CODE
  blk.alloc_frame_pos = 0;
#endif

  callstack = caml_callback(memprof_callback, Val_unit);
  blk.occurences = mt_generate_poisson(rest*lambda) + 1;

  if(caml_young_ptr - whsize < caml_young_trigger)
    caml_gc_dispatch();
  caml_young_ptr -= whsize;

  next_sample_young += whsize;
  update_limit();

  blk.block = Val_hp(caml_young_ptr);
  blk.callstack = callstack;
  push(&blk);

  CAMLreturn0;
}

value caml_memprof_track_alloc_shr(value block) {
  CAMLparam1(block);
  CAMLreturn (block);

  Assert(Is_in_heap(block));

  uint32_t occurences = mt_generate_poisson(lambda * Whsize_val(block));
  if(occurences > 0) {
    struct tracked_block blk;
    // We temporarily change the tag of the newly allocated  block to
    // Abstract_tag, because we may trigger the GC and we want to avoid
    // scanning this uninitialized block.
    tag_t oldtag = Tag_val(block);
    Tag_val(block) = Abstract_tag;

#ifdef NATIVE_CODE
    blk.alloc_frame_pos = 0;
#endif

    blk.occurences = occurences;
    blk.block = block;

    blk.callstack = caml_callback(memprof_callback, Val_unit);
    push(&blk);
    Tag_val(block) = oldtag;
  }

  CAMLreturn (block);
}

void caml_memprof_track_interned(header_t* block, header_t* blockend) {
  CAMLparam0();
  CAMLlocal1(root);
  struct tracked_block blk;

  value* sampled = NULL;
  double* rests = NULL;

  uintnat sz = 0;
  uintnat j = 0, i;
  header_t* p = block;

  /* We have to select the sampled blocks before inserting them,
     because insertion may trigger GC, and then blocks can escape
     from [block, blockend[ */
  while(1) {
    double next_sample = mt_generate_exponential();
    header_t *next_sample_p, *p0;
    if(next_sample >= blockend - p)
      break;

    p0 = p;
    next_sample_p = p + (uintnat)next_sample;
    Assert(next_sample_p <= blockend);
    while(p + Whsize_hp(p) <= next_sample_p)
      p += Whsize_hp(p);

    if(sz == 0) {
      sz = 32;
      sampled = (value*)malloc(sizeof(value) * sz);
      if(sampled == NULL)
        CAMLreturn0;
      rests = (double*)malloc(sizeof(double) * sz);
      if(rests == NULL) {
        free(sampled);
        CAMLreturn0;
      }
    } else if(j >= sz) {
      value* newsampled;
      double* newrests;
      sz *= 2;
      newsampled = (value*)realloc(sampled, sizeof(value) * sz);
      newrests = (double*)realloc(rests, sizeof(double) * sz);
      if(newsampled == NULL || newrests == NULL) {
        free(sampled);
        free(rests);
        CAMLreturn0;
      }
      sampled = newsampled;
      rests = newrests;
    }
    Assert (j < sz);

    sampled[j] = Val_hp(p);
    rests[j] = (p + Whsize_hp(p) - p0) - next_sample;
    Assert(rests[j] > 0);
    j++;

    p += Whsize_hp(p);
  }

  if(sz == 0)
    CAMLreturn0;

  CAMLxparamN(sampled, j);

#ifdef NATIVE_CODE
  blk.alloc_frame_pos = 0;
#endif

  /* Note : blk.callstack is not registered as a root, so we must take care
     not to call the GC before everything is pushed */
  blk.callstack = caml_callback_exn(memprof_callback, Val_unit);
  if (Is_exception_result(blk.callstack)) {
    free(sampled);
    free(rests);
    caml_raise(Extract_exception(blk.callstack));
  }

  for(i = 0; i < j; i++) {
    blk.occurences = mt_generate_poisson(rests[i]*lambda) + 1;
    blk.block = sampled[i];
    push(&blk);
  }

  free(sampled);
  free(rests);
  CAMLreturn0;
}

#ifdef NATIVE_CODE
double caml_memprof_call_gc_begin(void) {
  if(caml_young_ptr < caml_memprof_young_limit && lambda > 0) {
    double rest = caml_young_alloc_end - caml_young_ptr - next_sample_young;
    Assert(rest > 0);
    renew_sample();
    return rest;
  } else
    return 0;
}

void caml_memprof_call_gc_end(double rest) {
  uintnat pc = caml_last_return_address;
  char * sp = caml_bottom_of_stack;
  frame_descr * d;
  uintnat infoptr;
  uintnat i;
  uintnat tot_whsize;
  unsigned short num_blocks;
  unsigned short* block_sizes;
  CAMLparam0();
  CAMLlocal1(callstack);

  Assert(rest >= 0);

  if(lambda == 0)
    CAMLreturn0;
  d = caml_next_frame_descriptor(&pc, &sp);
  /* Should not happen, except on sparc */
  if(d == NULL || (d->frame_size & 2) != 2)
    CAMLreturn0;
  infoptr = ((uintnat) d +
             sizeof(char *) + sizeof(short) + sizeof(short) +
             sizeof(short) * d->num_live + sizeof(frame_descr *) - 1)
            & -sizeof(frame_descr *);

  num_blocks = *(unsigned short*)infoptr;
  block_sizes = ((unsigned short*)infoptr) + 1;

  tot_whsize = 0;
  for(i = 0; i < num_blocks; i++)
    tot_whsize += block_sizes[i];

  if(rest > 0) {
    double next_sample = tot_whsize - rest;
    uintnat offs = 0;
    struct tracked_block blk;

    callstack = caml_callback(memprof_callback, Val_unit);

    if(caml_young_ptr - tot_whsize < caml_young_trigger)
      caml_gc_dispatch();

    for(i = 0; i < num_blocks; i++) {
      next_sample -= block_sizes[i];
      offs += block_sizes[i];
      if(next_sample < 0) {
        blk.occurences = mt_generate_poisson(-next_sample*lambda) + 1;
        blk.alloc_frame_pos = i;

        /* It is important not to allocate anything or call the GC between
         * this point and the re-execution of the allocation code in assembly,
         * because blk.block would then be incorrect. */
        blk.block = Val_hp(caml_young_ptr - offs);
        blk.callstack = callstack;

        push(&blk);

        /* It is *not* garanteed that this block will actually be used,
         * because a signal can interrupt the allocation in the assembly code.
         * Thus, we put a special header on this block and we will lazily remove
         * this sample if this header is not overriden. */
        Hd_val(blk.block) = Make_header(0, 0, Caml_blue);

        next_sample = mt_generate_exponential();
      }
    }

    Assert(offs == tot_whsize);
  }

  /* We prevent the next allocation to be sampled, as it already had its
   * chance before the call. */
  next_sample_young += tot_whsize;
  update_limit();

  CAMLreturn0;
}
#endif

static void caml_memprof_clean() {
  uintnat i, j;

  for(i = j = old; i < tracked_blocks_end; i++) {
    value blk = tracked_blocks[i].block;
    Assert(Is_in_heap_or_young(blk));
    if(Is_young(blk)) {
      if(Hd_val(blk) == 0) {
        Assert(Is_in_heap(Field(blk, 0)));
        tracked_blocks[i].block = Field(blk, 0);
      } else if(Hd_val(blk) != Make_header(0, 0, Caml_blue)) {
        Assert(Is_white_val(blk) || Is_black_val(blk));
        Assert(Wosize_val(blk) != 0);
        Assert(Wosize_val(blk) <= Max_young_wosize);
      } else
        continue;
    }
    tracked_blocks[j] = tracked_blocks[i];
    j++;
  }
  tracked_blocks_end = j;

  shrink();
}

void caml_memprof_minor_gc_update(value* old_young_ptr) {
  uintnat i;

  for(i = old; i < tracked_blocks_end; i++) {
    value blk = tracked_blocks[i].block;
    Assert(Is_in_heap_or_young(blk));
    if(Is_young(blk)) {
      if(Hd_val(blk) == 0) {
        Assert(Is_in_heap(Field(blk, 0)));
        tracked_blocks[i].block = Field(blk, 0);
      } else {
        Assert(Hd_val(blk) == Make_header(0, 0, Caml_blue) ||
               Is_white_val(blk) || Is_black_val(blk));
        Assert(Hd_val(blk) == Make_header(0, 0, Caml_blue) ||
               Wosize_val(blk) != 0);
        Assert(Wosize_val(blk) <= Max_young_wosize);
        continue;
      }
    }
    tracked_blocks[old] = tracked_blocks[i];
    old++;
  }

  tracked_blocks_end = old;

  shrink();

  next_sample_young -= caml_young_alloc_end - old_young_ptr;
  update_limit();
}

void caml_memprof_do_young_roots(scanning_action f) {
  uintnat i;

  for(i = old; i < tracked_blocks_end; i++) {
    value* r = &tracked_blocks[i].callstack;
    f(*r, r);
  }
}

void caml_memprof_major_gc_update(void) {
  uintnat i, j;

  j = 0;
  for(i = 0; i < old; i++) {
    Assert(Is_in_heap(tracked_blocks[i].block));
    if(!Is_white_val(tracked_blocks[i].block))
      tracked_blocks[j++] = tracked_blocks[i];
  }
  old = j;
  for(; i < tracked_blocks_end; i++)
    if(Is_young(tracked_blocks[i].block) ||
       !Is_white_val(tracked_blocks[i].block))
      tracked_blocks[j++] = tracked_blocks[i];
  tracked_blocks_end = j;

  shrink();
}

void caml_memprof_do_strong_roots(scanning_action f) {
  uintnat i;

  for(i = 0; i < tracked_blocks_end; i++) {
    value* r = &tracked_blocks[i].callstack;
    f(*r, r);
  }
}

void caml_memprof_do_weak_roots(scanning_action f) {
  uintnat i;

  Assert(old == tracked_blocks_end);

  for(i = 0; i < tracked_blocks_end; i++) {
    value* r = &tracked_blocks[i].block;
    Assert(Is_in_heap(*r));
    f(*r, r);
  }
}

CAMLprim value caml_memprof_set(value v) {
  CAMLparam1(v);

  double l = Double_val(v);

  if(!(l >= 0) || l > 1)
    caml_invalid_argument("caml_memprof_set");

  lambda = l;
  lambda_rec = l == 0 ? INFINITY : 1/l;
  renew_sample();

  CAMLreturn(Val_unit);
}

CAMLprim value caml_memprof_clear(value unit) {
  CAMLparam0();

  tracked_blocks_end = 0;
  shrink();

  CAMLreturn(Val_unit);
}

CAMLprim value caml_memprof_dump_samples(value unit) {
  CAMLparam0();
  CAMLlocal3(res, sample, tmp);
  uintnat i;
  struct tracked_block* buffer;
  value* callstacks;

  caml_memprof_clean();

  if(tracked_blocks_end == 0)
    CAMLreturn(Atom(0)); // empty array

  if(tracked_blocks_end > Max_wosize)
    caml_failwith("Too many samples to dump.");

  buffer = (struct tracked_block*)
    malloc(tracked_blocks_end * sizeof(struct tracked_block));
  if(buffer == NULL)
    caml_raise_out_of_memory();
  memcpy(buffer, tracked_blocks,
         tracked_blocks_end * sizeof(struct tracked_block));

  callstacks = (value*)malloc(tracked_blocks_end * sizeof(value));
  if(callstacks == NULL) {
    free(buffer);
    caml_raise_out_of_memory();
  }

  for(i = 0; i < tracked_blocks_end; i++) {
    callstacks[i] = buffer[i].callstack;
    buffer[i].block = Val_long(Wosize_val(buffer[i].block));
  }

  Begin_roots_block(callstacks, tracked_blocks_end);
  res = caml_alloc(tracked_blocks_end, 0);

  for(i = 0; i < Wosize_val(res); i++) {
#ifdef NATIVE_CODE
    uintnat callstack_sz = Wosize_val(callstacks[i]);
#endif

    sample = caml_alloc(3, 0);
    Store_field(res, i, sample);
    Store_field(sample, 1, buffer[i].block);
    Store_field(sample, 2, Val_long(buffer[i].occurences));

#ifdef NATIVE_CODE
    if(buffer[i].alloc_frame_pos != 0 && callstack_sz > 1) {
      uintnat j;
      if (callstack_sz <= Max_young_wosize) {
        tmp = caml_alloc_small(callstack_sz, 0);
        for (j = 0; j < callstack_sz; j++)
          Field(tmp, j) = Field(callstacks[i], j);
        callstacks[i] = tmp;
      } else {
        tmp = caml_alloc_shr(callstack_sz, 0);
        for (j = 0; j < callstack_sz; j++)
          caml_initialize(&Field(tmp, j), Field(callstacks[i], j));
        callstacks[i] = tmp;
      }

      tmp = caml_alloc_small(2, 0);
      Field(tmp, 0) = Field(callstacks[i], 1);
      Field(tmp, 1) = Val_long(buffer[i].alloc_frame_pos);

      Store_field(callstacks[i], 1, tmp);
    }
#endif

    Store_field(sample, 0, callstacks[i]);
  }

  End_roots();

  free(buffer);
  free(callstacks);

  CAMLreturn(res);
}
