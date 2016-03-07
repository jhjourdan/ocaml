#include <math.h>
#include <string.h>
#include "caml/memprof.h"
#include "caml/backtrace_prim.h"
#include "caml/fail.h"
#include "caml/alloc.h"
#include "caml/callback.h"
#include "caml/weak.h"
#ifdef NATIVE_CODE
#include "stack.h"
#endif

static uint32_t mt_state[624];
static uint32_t mt_index;
static int mt_init = 0;

/* [lambda] is the mean number of samples for each allocated word (including
   block headers. */
static double lambda = 0;
static double lambda_rec = INFINITY;
static intnat callstack_size = 0;

static double next_sample_young;
value* caml_memprof_young_limit;

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

  caml_update_young_limit();
}

static void renew_sample(void) {
  next_sample_young =
    caml_young_alloc_end - caml_young_ptr + mt_generate_exponential();

  update_limit();
}

void caml_memprof_reinit(void) {
  int i;

  if(!mt_init) {
    mt_index = 624;
    mt_state[0] = 42;
    for(i = 1; i < 624; i++)
      mt_state[i] = 0x6c078965 * (mt_state[i-1] ^ (mt_state[i-1] >> 30)) + i;
  }

  renew_sample();
}

void set_lambda(double l) {
  lambda = l;
  lambda_rec = l == 0 ? INFINITY : 1/l;
  renew_sample();
}

CAMLprim value caml_memprof_set(value v) {
  CAMLparam1(v);
  double l = Double_val(Field(v, 0));
  intnat sz = Long_val(Field(v, 1));
  if(sz < 0 || !(l >= 0.) || l > 1.)
    caml_failwith("caml_memprof_set");
  set_lambda(l);
  callstack_size = sz;
  CAMLreturn(Val_unit);
}

CAMLprim value caml_memprof_get(value v) {
  CAMLparam1(v);
  CAMLlocal1(res);
  res = caml_alloc_small(2, 0);
  Field(res, 0) = caml_copy_double(lambda);
  Field(res, 1) = Val_long(callstack_size);
  CAMLreturn(res);
}

static value memprof_callback;

CAMLprim value caml_memprof_register_callback(value callback) {
  CAMLparam1(callback);
  memprof_callback = callback;
  caml_register_global_root(&memprof_callback);
  CAMLreturn(Val_unit);
}

static value capture_callstack() {
  CAMLparam0();
  CAMLlocal1(res);
  res = caml_get_current_callstack(Val_long(callstack_size));
  CAMLreturn(res);
}

static value do_callback(uintnat wosize, uintnat occurences, value callstack) {
  CAMLparam1(callstack);
  double old_lambda = lambda;
  set_lambda(0);
  value res = caml_callback3_exn(memprof_callback, Val_long(wosize),
                                 Val_long(occurences), callstack);
  set_lambda(old_lambda);
  CAMLreturn(res);
}

void caml_memprof_track_young(uintnat wosize) {
  CAMLparam0();
  CAMLlocal2(ephe, callstack);
  uintnat whsize = Whsize_wosize(wosize);
  double rest = caml_young_alloc_end - caml_young_ptr - next_sample_young;

  Assert(rest > 0);

  caml_young_ptr += whsize;

  long occurences = mt_generate_poisson(rest*lambda) + 1;
  callstack = capture_callstack();
  ephe = do_callback(wosize, occurences, callstack);
  if (Is_exception_result(ephe)) caml_raise(Extract_exception(ephe));

  if(caml_young_ptr - whsize < caml_young_trigger)
    caml_gc_dispatch();
  caml_young_ptr -= whsize;

  next_sample_young += whsize;
  update_limit();

  caml_ephe_set_key(ephe, Val_long(0), Val_hp(caml_young_ptr));

  CAMLreturn0;
}

value caml_memprof_track_alloc_shr(value block) {
  CAMLparam1(block);
  CAMLlocal2(ephe, callstack);

  Assert(Is_in_heap(block));

  uint32_t occurences = mt_generate_poisson(lambda * Whsize_val(block));
  if(occurences > 0) {
    // We temporarily change the tag of the newly allocated  block to
    // Abstract_tag, because we may trigger the GC and we want to avoid
    // scanning this uninitialized block.
    tag_t oldtag = Tag_val(block);
    Tag_val(block) = Abstract_tag;
    callstack = capture_callstack();
    ephe = do_callback(Wosize_val(block), occurences, callstack);
    if (Is_exception_result(ephe)) caml_raise(Extract_exception(ephe));
    Tag_val(block) = oldtag;

    caml_ephe_set_key(ephe, Val_long(0), block);
  }

  CAMLreturn (block);
}

void caml_memprof_track_interned(header_t* block, header_t* blockend) {
  CAMLparam0();
  CAMLlocal2(ephe, callstack);

  value* sampled = NULL;
  uintnat* occurences = NULL;
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
      occurences = (uintnat*)malloc(sizeof(uintnat) * sz);
      if(occurences == NULL) {
        free(sampled);
        CAMLreturn0;
      }
    } else if(j >= sz) {
      value* newsampled;
      uintnat* newoccurences;
      sz *= 2;
      newsampled = (value*)realloc(sampled, sizeof(value) * sz);
      if(newsampled == NULL) {
        free(sampled); free(occurences);
        CAMLreturn0;
      }
      sampled = newsampled;
      newoccurences = (uintnat*)realloc(occurences, sizeof(uintnat) * sz);
      if(newoccurences == NULL) {
        free(sampled); free(occurences);
        CAMLreturn0;
      }
      occurences = newoccurences;
    }
    Assert (j < sz);

    sampled[j] = Val_hp(p);
    double rest = (p + Whsize_hp(p) - p0) - next_sample;
    Assert(rest > 0);
    occurences[j] = mt_generate_poisson(rest*lambda) + 1;
    j++;

    p += Whsize_hp(p);
  }

  if(sz == 0)
    CAMLreturn0;

  CAMLxparamN(sampled, j);

  callstack = capture_callstack();
  for(i = 0; i < j; i++) {
    ephe = do_callback(Wosize_val(sampled[i]), occurences[i], callstack);
    if (Is_exception_result(ephe)) {
      free(sampled);
      free(occurences);
      caml_raise(Extract_exception(ephe));
    }
    caml_ephe_set_key(ephe, Val_long(0), sampled[i]);
  }

  free(sampled);
  free(occurences);
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
  uintnat i, j;
  uintnat tot_whsize;
  unsigned short num_blocks;
  unsigned short* block_sizes;
  CAMLparam0();
  CAMLlocal3(callstack, callstack_cur, tmp);

  Assert(rest >= 0);

  if(lambda == 0)
    CAMLreturn0;

  d = caml_next_frame_descriptor(&pc, &sp);
  Assert(d != NULL && (d->frame_size & 2) == 2);
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

    struct smp { uintnat offs, occurences, sz, alloc_frame_pos; } *samples, *p;
    uintnat sz = 2, n_samples = 0;
    samples = (struct smp*)malloc(sz*sizeof(struct smp));
    if(samples == NULL)
      goto abort;

    for(i = 0; i < num_blocks; i++) {
      next_sample -= block_sizes[i];
      offs += block_sizes[i];
      if(next_sample < 0) {
        if(n_samples >= sz) {
          sz *= 2;
          p = (struct smp*)realloc(samples, sz);
          if(p == NULL) { free(samples); goto abort; }
          samples = p;
        }

        samples[n_samples].offs = offs;
        samples[n_samples].occurences =
          mt_generate_poisson(-next_sample*lambda) + 1;
        samples[n_samples].sz = Wosize_whsize(block_sizes[i]);
        samples[n_samples].alloc_frame_pos = i;
        n_samples++;

        next_sample = mt_generate_exponential();
      }
    }

    Assert(offs == tot_whsize);

    CAMLlocalN(ephes, n_samples);
    callstack = capture_callstack();
    for(i = 0; i < n_samples; i++) {
      if(samples[i].alloc_frame_pos == 0)
        callstack_cur = callstack;
      else {
        callstack_cur = caml_alloc(Wosize_val(callstack), 0);
        for(j = 1; j < Wosize_val(callstack); j++)
          Store_field(callstack_cur, j, Field(callstack, j));
        tmp = caml_alloc_small(2, 0);
        Field(tmp, 0) = Field(callstack, 0);
        Field(tmp, 1) = Val_long(samples[i].alloc_frame_pos);
        Store_field(callstack_cur, 0, tmp);
      }

      ephes[i] = do_callback(samples[i].sz, samples[i].occurences, callstack_cur);
      if(Is_exception_result(ephes[i])) {
        free(samples);
        caml_raise(Extract_exception(ephes[i]));
      }
    }

    if(caml_young_ptr - tot_whsize < caml_young_trigger)
      caml_gc_dispatch();

    for(i = 0; i < n_samples; i++) {
      value v = Val_hp(caml_young_ptr - samples[i].offs);

      /* It is important not to allocate anything or call the GC between
       * this point and the re-execution of the allocation code in assembly,
       * because blk.block would then be incorrect. */
      caml_ephe_set_key(ephes[i], Val_long(0), v);

      /* It is *not* garanteed that this block will actually be used,
       * because a signal can interrupt the allocation in the assembly code.
       * Thus, we put an abstract header on this unitialized block. */
      Hd_val(v) = Make_header(samples[i].sz, Abstract_tag, Caml_black);
    }

    free(samples);
  }

  /* We prevent the next allocation to be sampled, as it already had its
   * chance before the call. */
 abort:
  next_sample_young += tot_whsize;
  update_limit();

  CAMLreturn0;
}
#endif

void caml_memprof_minor_gc_update(value* old_young_ptr) {
  next_sample_young -= caml_young_alloc_end - old_young_ptr;
  update_limit();
}
