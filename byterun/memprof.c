#include <math.h>
#include "memprof.h"
#include "backtrace.h"
#include "signals.h"
#include "stacks.h"

static uint32 mt_state[624];
static uint32 mt_index;

/* [lambda] is the sampling rate. */
static double lambda = 0.01;

static double next_sample_pos;
char* caml_memprof_young_limit;

struct caml_memprof_tracked_block* caml_memprof_tracked_blocks = NULL;
uintnat caml_memprof_tracked_blocks_end = 0;
static uintnat size = 0, old = 0;

static double mt_generate_uniform(void) {
  int i;
  uint32 y;

  /* Mersenne twister PRNG */
  if (mt_index == 0)
    for(i = 0; i < 624; i++) {
      y = (mt_state[i] & 0x80000000) + (mt_state[(i+1)%624] & 0x7fffffff);
      mt_state[i] = mt_state[(i+397)%624] ^ (y >> 1);
      if (y&1)
        mt_state[i] ^= 0x9908b0df;
    }

  y = mt_state[mt_index];
  y = y ^ (y >> 11);
  y = y ^ ((y << 7) & 0x9d2c5680);
  y = y ^ ((y << 15) & 0xefc60000);
  y = y ^ (y >> 18);

  mt_index = (mt_index + 1) % 624;
  return y*2.3283064365386962890625e-10 + /* 2^-32 */
          1.16415321826934814453125e-10; /* 2^-33 */
}

static double mt_generate_exponential(double lambda) {
  Assert(lambda >= 0 && !isinf(lambda));

  if(lambda == 0)
    return INFINITY;

  return -logf(mt_generate_uniform())/lambda;
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
    k = logf(c) - lambda - logf(beta);
    while(1) {
      double u, x, n, v, y, y2;
      u = mt_generate_uniform();
      x = (alpha - logf((1.0 - u)/u))/beta;
      n = floorf(x + 0.5);
      if(n < 0.)
        continue;

      v = mt_generate_uniform();
      y = alpha - beta*x;
      y2 = 1. + expf(y);
      if(y + logf(v/(y2*y2)) < k + n*logf(lambda) - lgammaf(n+1))
        return (n > ((uint32)-2)) ? ((uint32)-2) : n;
    }
  }
}

static void update_limit(void) {
  if(lambda == 0)
    caml_memprof_young_limit = caml_young_start;
  else
    caml_memprof_young_limit =
      caml_young_end - Bsize_wsize((uintnat)next_sample_pos);

  intnat oldpending = caml_signals_are_pending;
  if(caml_young_limit < caml_memprof_young_limit)
    caml_young_limit = caml_memprof_young_limit;
  if(caml_signals_are_pending && !oldpending)
    caml_young_limit = caml_young_end;
}

static void renew_sample(void) {
  next_sample_pos =
    Wsize_bsize(caml_young_end-caml_young_ptr) + mt_generate_exponential(lambda);

  update_limit();

  Assert(caml_memprof_young_limit <= caml_young_ptr ||
         caml_memprof_young_limit == caml_young_start);
}

void caml_memprof_reinit(void) {
  int i;

  if(caml_memprof_tracked_blocks == NULL) {
    mt_index = 0;
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

  renew_sample();
}

static void caml_memprof_tracked_blocks_add(struct caml_memprof_tracked_block* blk,
                                            uintnat num) {
  uintnat newsize = size;
  uintnat i;

  while(caml_memprof_tracked_blocks_end+num > newsize) newsize *= 2;

  if(size != newsize) {
    caml_memprof_tracked_blocks =
      (struct caml_memprof_tracked_block*)realloc(caml_memprof_tracked_blocks,
         newsize * sizeof(struct caml_memprof_tracked_block));
    if(caml_memprof_tracked_blocks == NULL)
      caml_fatal_error("out of memory");
    size = newsize;
  }

  Assert(caml_memprof_tracked_blocks_end + num < size);

  for(i = 0; i < num; i++) {
    Assert(Is_in_heap_or_young(blk[i].block));
    caml_memprof_tracked_blocks[caml_memprof_tracked_blocks_end] = blk[i];
    caml_memprof_tracked_blocks_end++;
  }
}

static void caml_memprof_tracked_blocks_shrink(void) {
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

void caml_memprof_track_one(value block, uintnat wosize) {
  double rest;

  Assert(Is_in_heap_or_young(block));

  if(Is_in_heap(block)) {
    next_sample_pos -= Whsize_wosize(wosize);
    update_limit();
  }

  rest = Wsize_bsize(caml_young_end-caml_young_ptr) - next_sample_pos;
  if(rest > 0) {
    struct caml_memprof_tracked_block blk;

    renew_sample();

#ifdef NATIVE_CODE
    {
      uintnat pc = caml_last_return_address;
      char * sp = caml_bottom_of_stack;
      frame_descr * d;
      d = caml_next_frame_descriptor(&pc, &sp);
      if(d == NULL || (d -> frame_size & 3) != 1)
        blk.loc1 = blk.loc2 = 0;
      else {
        uintnat infoptr =
          ((uintnat) d +
           sizeof(char *) + sizeof(short) + sizeof(short) +
           sizeof(short) * d->num_live + sizeof(frame_descr *) - 1)
          & -sizeof(frame_descr *);
        blk.loc1 = ((uint32 *)infoptr)[0];
        blk.loc2 = ((uint32 *)infoptr)[1];
      }
    }
#else
    {
      value * sp = caml_extern_sp;
      value * trapsp = caml_trapsp;
      blk.loc = caml_next_frame_pointer(&sp, &trapsp);
    }
#endif

    blk.occurences = mt_generate_poisson(rest*lambda) + 1;
    blk.block = block;
    caml_memprof_tracked_blocks_add(&blk, 1);
  }
}

#ifdef NATIVE_CODE
double caml_memprof_call_gc_begin_hook(void) {
  if(caml_young_ptr < caml_memprof_young_limit && lambda > 0) {
    double rest = Wsize_bsize(caml_young_end-caml_young_ptr) - next_sample_pos;
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

void caml_memprof_call_gc_end_hook(double rest) {
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
      caml_fatal_error("out of memory");

    for(i = 0; i < num_blocks; i++) {
      next_sample -= dbg[i].whsize;
      offs += dbg[i].whsize;
      if(next_sample < 0) {
        blk[blk_pos].occurences = mt_generate_poisson(-next_sample*lambda) + 1;
        next_sample = mt_generate_exponential(lambda);

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

    caml_memprof_tracked_blocks_add(blk, blk_pos);

    free(blk);
  }

  /* We prevent the next allocation to be sampled, as it already had its
   * chance before the call. */
  next_sample_pos += tot_whsize;
  update_limit();
}
#endif

void caml_memprof_minor_gc() {
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

  caml_memprof_tracked_blocks_shrink();

  renew_sample();
}

void caml_memprof_clean() {
  uintnat i;

  for(i = old; i < caml_memprof_tracked_blocks_end; i++) {
    value blk = caml_memprof_tracked_blocks[i].block;
    Assert(Is_in_heap_or_young(blk));
    if(Is_young(blk)) {
      if(Hd_val(blk) == 0) {
        Assert(Is_in_heap(Field(blk, 0)));
        caml_memprof_tracked_blocks[i].block = Field(blk, 0);
      } else if(Hd_val(blk) == Make_header(0, 0, Caml_blue))
        continue;
      else {
        Assert(Is_white_val(blk) || Is_black_val(blk));
        Assert(Wosize_val(blk) != 0);
        Assert(Wosize_val(blk) <= Max_young_wosize);
      }
    }
    caml_memprof_tracked_blocks[old] = caml_memprof_tracked_blocks[i];
    old++;
  }
  caml_memprof_tracked_blocks_end = old;

  caml_memprof_tracked_blocks_shrink();
}

void caml_memprof_major_gc(void) {
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

  caml_memprof_tracked_blocks_shrink();
}

void caml_memprof_do_weak_roots(scanning_action f) {
  uintnat i;

  Assert(old == caml_memprof_tracked_blocks_end);

  for(i = 0; i < caml_memprof_tracked_blocks_end; i++) {
    Assert(Is_in_heap(caml_memprof_tracked_blocks[i].block));
    f(caml_memprof_tracked_blocks[i].block,
      &caml_memprof_tracked_blocks[i].block);
  }
}

double caml_memprof_get_lambda(void) {
  return lambda;
}


void caml_memprof_set_lambda(double l) {
  lambda = l;
  renew_sample();
}
