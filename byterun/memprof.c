#include <math.h>
#include "memprof.h"
#include "backtrace.h"
#include "signals.h"
#include "stacks.h"
#include "fail.h"
#include "memory.h"
#include "alloc.h"
#include "memory.h"

static uint32 mt_state[624];
static uint32 mt_index;

/* [lambda] is the mean number of samples for each allocated word (including
   block headers. */
static double lambda = 0.01;
static double lambda_rec = 100;

static intnat dumpped_backtrace_size = 3;

static double next_sample_young;
char* caml_memprof_young_limit;

struct caml_memprof_tracked_block* caml_memprof_tracked_blocks = NULL;
uintnat caml_memprof_tracked_blocks_end = 0;
static uintnat size = 0, old = 0;

static float log_tbl[64];

static double fast_log(float v) {
  int32 iv = *(int32*)&v;
  int32 exp = ((iv >> 23) & 255) - 127;
  int32 idx = iv & 0x7E0000;
  double corr = 1. * (0x1FFFF & iv) / (idx + 0x800000);

  return ((double)exp + log_tbl[idx >> 17]) * 0.69314718055994530941 + corr;
}

static double mt_generate_uniform(void) {
  int i;
  uint32 y;

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

  double res = -fast_log(mt_generate_uniform()) * lambda_rec;
  if(res < 0) return 0;
  return res;
}

/* Max returned value : 2^32-2 */
static uint32 mt_generate_poisson(double lambda) {
  Assert(lambda >= 0 && !isinf(lambda));

  if(lambda == 0)
    return 0;

  if(lambda < 20) {
    double p;
    uint32 k;
    k = 0;
    p = expf(lambda);
    do {
      k++;
      p *= mt_generate_uniform();
    } while(p > 1);
    return k-1;
  } else {
    double c, beta, alpha, k;
    c = 0.767 - 3.36/lambda;
    beta = M_PI/sqrt(3.*lambda);
    alpha = beta*lambda;
    k = fast_log(c) - lambda - fast_log(beta);
    while(1) {
      double u, x, n, v, y, y2;
      u = mt_generate_uniform();
      x = (alpha - fast_log((1.0 - u)/u))/beta;
      n = floor(x + 0.5);
      if(n < 0.)
        continue;

      v = mt_generate_uniform();
      y = alpha - beta*x;
      y2 = 1. + expf(y);
      if(y + fast_log(v/(y2*y2)) < k + n*fast_log(lambda) - lgammaf(n+1))
        return (n > ((uint32)-2)) ? ((uint32)-2) : n;
    }
  }
}

static void update_limit(void) {
  if(next_sample_young >
     Wsize_bsize(caml_young_end - caml_young_start) +
       Whsize_wosize(Max_young_wosize))
    caml_memprof_young_limit = caml_young_start - Max_young_wosize;
  else
    caml_memprof_young_limit =
      (char*)caml_young_end - Bsize_wsize((uintnat)next_sample_young);

  intnat oldpending = caml_signals_are_pending;

  if(caml_young_limit != caml_young_end)
    caml_young_limit =
      caml_memprof_young_limit < caml_young_start ?
      caml_young_start : caml_memprof_young_limit;
  if(caml_signals_are_pending && !oldpending)
    caml_young_limit = caml_young_end;
}

static void renew_sample(void) {
  next_sample_young =
    Wsize_bsize(caml_young_end-caml_young_ptr) + mt_generate_exponential();

  update_limit();

  Assert(caml_memprof_young_limit <= caml_young_ptr ||
         caml_memprof_young_limit == caml_young_start);
}

void caml_memprof_reinit(void) {
  int i;

  if(caml_memprof_tracked_blocks == NULL) {
    mt_index = 624;
    mt_state[0] = 42;
    for(i = 1; i < 624; i++)
      mt_state[i] = 0x6c078965 * (mt_state[i-1] ^ (mt_state[i-1] >> 30)) + i;

    caml_memprof_tracked_blocks =
      (struct caml_memprof_tracked_block*)
      malloc(1024 * sizeof(struct caml_memprof_tracked_block));
    if(caml_memprof_tracked_blocks == NULL)
      caml_fatal_error("out of memory");
    size = 1024;
  }

  for(i = 0; i < 64; ++i)
    log_tbl[i] =
      (log(1. + (2*i+1) / 128.) - 1. / (2*i + 128)) / log(2.);

  renew_sample();
}

static void push_blks(struct caml_memprof_tracked_block* blk, uintnat num) {
  uintnat newsize = size;
  uintnat i;

  while(caml_memprof_tracked_blocks_end + num > newsize) newsize *= 2;

  if(size != newsize) {
    struct caml_memprof_tracked_block* p =
      (struct caml_memprof_tracked_block*)realloc(caml_memprof_tracked_blocks,
         newsize * sizeof(struct caml_memprof_tracked_block));

    if(p == NULL)
      return;

    caml_memprof_tracked_blocks = p;

    size = newsize;
  }

  Assert(caml_memprof_tracked_blocks_end + num <= size);

  for(i = 0; i < num; i++) {
    Assert(Is_in_heap_or_young(blk[i].block));
    caml_memprof_tracked_blocks[caml_memprof_tracked_blocks_end] = blk[i];
    caml_memprof_tracked_blocks_end++;
    caml_remove_generational_global_root(&blk[i].callstack);
  }
}

static void shrink(void) {
  uintnat newsize = size;
  while(newsize >= 2048 && caml_memprof_tracked_blocks_end * 8 < newsize)
    newsize /= 2;
  if(newsize != size) {
    struct caml_memprof_tracked_block* p =
      (struct caml_memprof_tracked_block*)realloc(caml_memprof_tracked_blocks,
         newsize * sizeof(struct caml_memprof_tracked_block));
    if(p != NULL) {
      caml_memprof_tracked_blocks = p;
      size = newsize;
    }
  }
}

static void fill_blk(struct caml_memprof_tracked_block *blk) {
  #ifdef NATIVE_CODE
  blk -> loc1 = blk -> loc2 = 0;
  #endif

  blk -> callstack =
    //caml_alloc((mlsize_t) 12, Abstract_tag);
    caml_get_current_callstack(Val_long(dumpped_backtrace_size));
  caml_register_generational_global_root(&(blk -> callstack));
}

void caml_memprof_track_young(uintnat wosize) {
  double rest = Wsize_bsize(caml_young_end-caml_young_ptr) - next_sample_young;
  struct caml_memprof_tracked_block blk;

  Assert(rest > 0);

  caml_young_ptr += Bhsize_wosize(wosize);

  renew_sample();

  fill_blk(&blk);
  blk.occurences = mt_generate_poisson(rest*lambda) + 1;

  if(caml_young_ptr - Bhsize_wosize(wosize) < caml_young_start)
    caml_minor_collection();
  caml_young_ptr -= Bhsize_wosize(wosize);

  next_sample_young += Whsize_wosize(wosize);
  update_limit();

  blk.block = Val_hp(caml_young_ptr);
  push_blks(&blk, 1);
}

value caml_memprof_track_alloc_shr(value block) {
  CAMLparam1(block);

  Assert(Is_in_heap(block));

  uint32 occurences = mt_generate_poisson(lambda * Whsize_val(block));
  if(occurences > 0) {
    struct caml_memprof_tracked_block blk;
    tag_t oldtag = Tag_val(block);
    Tag_val(block) = Abstract_tag;

    fill_blk(&blk);
    blk.occurences = occurences;
    blk.block = block;

    push_blks(&blk, 1);
    Tag_val(block) = oldtag;
  }

  CAMLreturn (block);
}

void caml_memprof_track_interned(header_t* block, header_t* blockend) {
  CAMLparam0();
  CAMLlocal1(root);

  value* sampled = NULL;
  double* rests;

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

  CAMLxparamN(sampled, j);

  for(i = 0; i < j; i++) {
    struct caml_memprof_tracked_block blk;

    fill_blk(&blk);
    blk.occurences = mt_generate_poisson(rests[i]*lambda) + 1;
    blk.block = sampled[i];
    push_blks(&blk, 1);
  }

  if(sz > 0) {
    free(sampled);
    free(rests);
  }
  CAMLreturn0;
}

#ifdef NATIVE_CODE
double caml_memprof_call_gc_begin(void) {
  if(caml_young_ptr < caml_memprof_young_limit && lambda > 0) {
    double rest = Wsize_bsize(caml_young_end-caml_young_ptr) - next_sample_young;
    Assert(rest > 0);
    renew_sample();
    return rest;
  } else
    return 0;
}

struct caml_memprof_alloc_debuginfo {
  uintnat whsize;
  uint32 loc1, loc2;
};

void caml_memprof_call_gc_end(double rest) {
  uintnat pc = caml_last_return_address;
  char * sp = caml_bottom_of_stack;
  frame_descr * d;
  uintnat infoptr;
  uintnat i;
  uintnat tot_whsize;
  uintnat num_blocks;
  struct caml_memprof_alloc_debuginfo* dbg;

  Assert(rest >= 0);

  if(lambda == 0)
    return;

  d = caml_next_frame_descriptor(&pc, &sp);
  if(d == NULL || (d -> frame_size & 3) != 2)
    caml_fatal_error("unable to find alloc frame descriptor");
  infoptr = ((uintnat) d +
             sizeof(char *) + sizeof(short) + sizeof(short) +
             sizeof(short) * d->num_live + sizeof(frame_descr *) - 1)
            & -sizeof(frame_descr *);

  num_blocks = *(uintnat*)infoptr;
  dbg = (struct caml_memprof_alloc_debuginfo*)(infoptr+sizeof(uintnat));

  tot_whsize = 0;
  for(i = 0; i < num_blocks; i++)
    tot_whsize += dbg[i].whsize;

  if(rest > 0) {
    double next_sample = tot_whsize - rest;
    uintnat blk_pos = 0;
    struct caml_memprof_tracked_block* blk =
      malloc(num_blocks * sizeof(struct caml_memprof_tracked_block));
    uintnat offs = 0;
    if(blk == NULL)
      return;

    for(i = 0; i < num_blocks; i++) {
      next_sample -= dbg[i].whsize;
      offs += dbg[i].whsize;
      if(next_sample < 0) {
        fill_blk(&blk[blk_pos]);

        blk[blk_pos].occurences = mt_generate_poisson(-next_sample*lambda) + 1;
        next_sample = mt_generate_exponential();

        blk[blk_pos].block = Val_long(offs);

        blk[blk_pos].loc1 = dbg[i].loc1;
        blk[blk_pos].loc2 = dbg[i].loc2;

        blk_pos++;
      }
    }

    Assert(offs == tot_whsize);

    if(caml_young_ptr - Bsize_wsize(tot_whsize) < caml_young_start)
      caml_minor_collection();

    for(i = 0; i < blk_pos; i++) {
      blk[i].block = Val_hp(caml_young_ptr -
                            Bsize_wsize(Unsigned_long_val(blk[i].block)));

      /* It is *not* garanteed that this block will actually be used,
       * because a signal can interrupt the allocation in the assembly code.
       * Thus, we put a special header on this block and we will lazily remove
       * this sample it this header is not overriden. */
      Hd_val(blk[i].block) = Make_header(0, 0, Caml_blue);
    }

    push_blks(blk, blk_pos);

    free(blk);
  }

  /* We prevent the next allocation to be sampled, as it already had its
   * chance before the call. */
  next_sample_young += tot_whsize;
  update_limit();
}
#endif

void caml_memprof_clean() {
  uintnat i, j;

  for(i = j = old; i < caml_memprof_tracked_blocks_end; i++) {
    value blk = caml_memprof_tracked_blocks[i].block;
    Assert(Is_in_heap_or_young(blk));
    if(Is_young(blk)) {
      if(Hd_val(blk) == 0) {
        Assert(Is_in_heap(Field(blk, 0)));
        caml_memprof_tracked_blocks[i].block = Field(blk, 0);
      } else if(Hd_val(blk) != Make_header(0, 0, Caml_blue)) {
        Assert(Is_white_val(blk) || Is_black_val(blk));
        Assert(Wosize_val(blk) != 0);
        Assert(Wosize_val(blk) <= Max_young_wosize);
      } else
        continue;
    }
    caml_memprof_tracked_blocks[j] = caml_memprof_tracked_blocks[i];
    j++;
  }
  caml_memprof_tracked_blocks_end = j;

  shrink();
}

void caml_memprof_minor_gc_update(char* old_young_ptr) {
  uintnat i;

  for(i = old; i < caml_memprof_tracked_blocks_end; i++) {
    value blk = caml_memprof_tracked_blocks[i].block;
    Assert(Is_in_heap_or_young(blk));
    if(Is_young(blk)) {
      if(Hd_val(blk) == 0) {
        Assert(Is_in_heap(Field(blk, 0)));
        caml_memprof_tracked_blocks[i].block = Field(blk, 0);
      } else {
        Assert(Hd_val(blk) == Make_header(0, 0, Caml_blue) ||
               Is_white_val(blk) || Is_black_val(blk));
        Assert(Hd_val(blk) == Make_header(0, 0, Caml_blue) ||
               Wosize_val(blk) != 0);
        Assert(Wosize_val(blk) <= Max_young_wosize);
        continue;
      }
    }
    caml_memprof_tracked_blocks[old] = caml_memprof_tracked_blocks[i];
    old++;
  }

  caml_memprof_tracked_blocks_end = old;

  shrink();

  next_sample_young -= Wsize_bsize(caml_young_end-old_young_ptr);
  update_limit();
}

void caml_memprof_do_young_roots(scanning_action f) {
  uintnat i;

  for(i = old; i < caml_memprof_tracked_blocks_end; i++) {
    value* r = &caml_memprof_tracked_blocks[i].callstack;
    f(*r, r);
  }
}

void caml_memprof_major_gc_update(void) {
  uintnat i, j;

  Assert(old == caml_memprof_tracked_blocks_end);

  j = 0;
  for(i = 0; i < caml_memprof_tracked_blocks_end; i++) {
    Assert(Is_in_heap(caml_memprof_tracked_blocks[i].block));
    if(!Is_white_val(caml_memprof_tracked_blocks[i].block))
      caml_memprof_tracked_blocks[j++] =
        caml_memprof_tracked_blocks[i];
  }

  old = caml_memprof_tracked_blocks_end = j;

  shrink();
}

void caml_memprof_do_strong_roots(scanning_action f) {
  uintnat i;

  Assert(old == caml_memprof_tracked_blocks_end);

  for(i = 0; i < caml_memprof_tracked_blocks_end; i++) {
    value* r = &caml_memprof_tracked_blocks[i].callstack;
    f(*r, r);
  }
}

void caml_memprof_do_weak_roots(scanning_action f) {
  uintnat i;

  Assert(old == caml_memprof_tracked_blocks_end);

  for(i = 0; i < caml_memprof_tracked_blocks_end; i++) {
    value* r = &caml_memprof_tracked_blocks[i].block;
    Assert(Is_in_heap(*r));
    f(*r, r);
  }
}

CAMLprim value caml_memprof_get(value unit) {
  CAMLparam0();
  CAMLlocal1(res);

  res = caml_alloc(2, 0);
  Store_field(res, 0, caml_copy_double(lambda));
  Store_field(res, 1, Val_long(dumpped_backtrace_size));

  CAMLreturn(res);
}


CAMLprim value caml_memprof_set(value v) {
  CAMLparam1(v);

  double l = Double_val(Field(v, 0));
  intnat sz = Long_val(Field(v, 1));

  if(!(l >= 0) || l > 1 || sz < 0)
    caml_invalid_argument("caml_memprof_set");

  lambda = l;
  lambda_rec = l == 0 ? INFINITY : 1/l;
  renew_sample();

  dumpped_backtrace_size = sz;

  CAMLreturn(Val_unit);
}

CAMLprim value caml_memprof_clear(value unit) {
  CAMLparam0();

  caml_memprof_tracked_blocks_end = 0;
  shrink();

  CAMLreturn(Val_unit);
}

CAMLprim value caml_memprof_estimate_consumed_memory(value unit) {
  CAMLparam0();
  uintnat i;
  intnat tot = 0;

  caml_memprof_clean();
  for(i = 0; i < caml_memprof_tracked_blocks_end; i++)
    tot += caml_memprof_tracked_blocks[i].occurences;

  tot = lrint(tot * lambda_rec);

  CAMLreturn(Val_long(tot));
}
