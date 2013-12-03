#include <math.h>
#include <string.h>
#include "memprof.h"
#include "backtrace.h"
#include "signals.h"
#include "stacks.h"
#include "fail.h"
#include "memory.h"
#include "alloc.h"
#include "hash.h"

static uint32 mt_state[624];
static uint32 mt_index;

/* [lambda] is the mean number of samples for each allocated word (including
   block headers. */
static double lambda = 0;
static double lambda_rec = INFINITY;

static intnat dumpped_backtrace_size = 3;

static double next_sample_young;
char* caml_memprof_young_limit;

struct tracked_block {
  value block;
#ifdef NATIVE_CODE
  uint32 alloc_frame_pos;
#endif
  uint32 occurences;
  value callstack;
};

struct tracked_block* tracked_blocks = NULL;
uintnat tracked_blocks_end = 0;
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

/* Max returned value : 2^30-2 */
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
        return n > ((1<<30)-2) ? ((1<<30)-2) : n;
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

  for(i = 0; i < 64; ++i)
    log_tbl[i] =
      (log(1. + (2*i+1) / 128.) - 1. / (2*i + 128)) / log(2.);

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

void caml_memprof_track_young(uintnat wosize) {
  CAMLparam0();
  CAMLlocal1(callstack);
  double rest = Wsize_bsize(caml_young_end-caml_young_ptr) - next_sample_young;
  struct tracked_block blk;

  Assert(rest > 0);

  caml_young_ptr += Bhsize_wosize(wosize);

  renew_sample();

#ifdef NATIVE_CODE
  blk.alloc_frame_pos = 0;
#endif

  callstack = caml_get_current_callstack(Val_long(dumpped_backtrace_size));
  blk.occurences = mt_generate_poisson(rest*lambda) + 1;

  if(caml_young_ptr - Bhsize_wosize(wosize) < caml_young_start)
    caml_minor_collection();
  caml_young_ptr -= Bhsize_wosize(wosize);

  next_sample_young += Whsize_wosize(wosize);
  update_limit();

  blk.block = Val_hp(caml_young_ptr);
  blk.callstack = callstack;
  push(&blk);

  CAMLreturn0;
}

value caml_memprof_track_alloc_shr(value block) {
  CAMLparam1(block);

  Assert(Is_in_heap(block));

  uint32 occurences = mt_generate_poisson(lambda * Whsize_val(block));
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

    blk.callstack = caml_get_current_callstack(Val_long(dumpped_backtrace_size));
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

#ifdef NATIVE_CODE
  blk.alloc_frame_pos = 0;
#endif

  /* Note : blk.callstack is not registered as a root, so we must take care
     not to call the GC before everything is pushed */
  blk.callstack = caml_get_current_callstack(Val_long(dumpped_backtrace_size));

  for(i = 0; i < j; i++) {
    blk.occurences = mt_generate_poisson(rests[i]*lambda) + 1;
    blk.block = sampled[i];
    push(&blk);
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

    callstack = caml_get_current_callstack(Val_long(dumpped_backtrace_size));

    if(caml_young_ptr - Bsize_wsize(tot_whsize) < caml_young_start)
      caml_minor_collection();

    for(i = 0; i < num_blocks; i++) {
      next_sample -= block_sizes[i];
      offs += block_sizes[i];
      if(next_sample < 0) {
        blk.occurences = mt_generate_poisson(-next_sample*lambda) + 1;
        blk.alloc_frame_pos = i;

        /* It is important not to allocate anything or call the GC between
         * this point and the re-execution of the allocation code in assembly,
         * because blk.block would then be incorrect. */
        blk.block = Val_hp(caml_young_ptr - Bsize_wsize(offs));
        blk.callstack = callstack;

        push(&blk);

        /* It is *not* garanteed that this block will actually be used,
         * because a signal can interrupt the allocation in the assembly code.
         * Thus, we put a special header on this block and we will lazily remove
         * this sample it this header is not overriden. */
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

void caml_memprof_clean() {
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

void caml_memprof_minor_gc_update(char* old_young_ptr) {
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

  next_sample_young -= Wsize_bsize(caml_young_end-old_young_ptr);
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

  Assert(old == tracked_blocks_end);

  j = 0;
  for(i = 0; i < tracked_blocks_end; i++) {
    Assert(Is_in_heap(tracked_blocks[i].block));
    if(!Is_white_val(tracked_blocks[i].block))
      tracked_blocks[j++] = tracked_blocks[i];
  }

  old = tracked_blocks_end = j;

  shrink();
}

void caml_memprof_do_strong_roots(scanning_action f) {
  uintnat i;

  Assert(old == tracked_blocks_end);

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

  tracked_blocks_end = 0;
  shrink();

  CAMLreturn(Val_unit);
}

static value value_of_locinfo(const struct loc_info *li, int32 h) {
  CAMLparam0();
  CAMLlocal1(res);

  res = caml_alloc(5, 0);

  Store_field(res, 0, Val_long(h & 0x3FFFFFFF));

  if(li->loc_valid) {
    Store_field(res, 1, caml_copy_string(li->loc_filename));
    Store_field(res, 2, Val_long(li->loc_lnum));
    Store_field(res, 3, Val_long(li->loc_startchr));
    Store_field(res, 4, Val_long(li->loc_endchr));
  } else {
    Store_field(res, 1, caml_copy_string(""));
    Store_field(res, 2, Val_long(0));
    Store_field(res, 3, Val_long(0));
    Store_field(res, 4, Val_long(0));
  }

  CAMLreturn(res);
}

struct loc_htbl {
  struct loc_entry {
    intnat isample, jloc;
    uint32 h;
  } *htbl;
  uintnat lsize, slots;
};

static void loc_htbl_resize(struct loc_htbl* loc_htbl) {
  uintnat i, new_lsize, mask;
  struct loc_entry *new_htbl;

  if(loc_htbl->htbl != NULL &&
     loc_htbl->slots <= (1LL << (loc_htbl->lsize - 1)))
    return;

  new_lsize = loc_htbl->htbl == NULL ? 6 : loc_htbl->lsize + 1;
  mask = (1 << new_lsize) - 1;
  new_htbl =
    (struct loc_entry*)malloc(sizeof(struct loc_entry) << new_lsize);

  if(new_htbl == NULL) {
    free(loc_htbl->htbl);
    loc_htbl->htbl = NULL;
    return;
  }

  for(i = 0; (i >> new_lsize) == 0; i++)
    new_htbl[i].isample = -1;

  if(loc_htbl->htbl != NULL) {
    for(i = 0; (i >> loc_htbl->lsize) == 0; i++) {
      uintnat j = loc_htbl->htbl[i].h & mask;
      while(new_htbl[j].isample != -1) j = (j+1) & mask;
      new_htbl[j] = loc_htbl->htbl[i];
    }
    free(loc_htbl->htbl);
  }

  loc_htbl->htbl = new_htbl;
  loc_htbl->lsize = new_lsize;
}

CAMLprim value caml_memprof_dump_samples(value unit) {
  CAMLparam0();
  CAMLlocal4(res, sample, trace, events);
  uintnat i, j;
  struct tracked_block* buffer;
  struct loc_htbl loc_htbl;
  value* callstacks;

  caml_memprof_clean();

  if(tracked_blocks_end == 0)
    CAMLreturn(Atom(0)); // empty array

  if(tracked_blocks_end > Max_wosize)
    caml_failwith("Too many samples to dump.");

#ifndef NATIVE_CODE
  events = caml_read_debug_info();
  if(events == Val_false)
    caml_failwith(caml_read_debug_info_error);
#endif

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

  loc_htbl.htbl = NULL;
  loc_htbl.slots = 0;
  loc_htbl_resize(&loc_htbl);
  if(loc_htbl.htbl == NULL) goto oom;

  for(i = 0; i < Wosize_val(res); i++) {
    struct loc_info li;

    sample = caml_alloc(3, 0);
    Store_field(res, i, sample);
    Store_field(sample, 1, buffer[i].block);
    Store_field(sample, 2, Val_long(buffer[i].occurences));

    trace = caml_alloc(Wosize_val(callstacks[i]), 0);
    Store_field(sample, 0, trace);

    for(j = 0; j < Wosize_val(callstacks[i]); j++) {
#ifdef NATIVE_CODE
      frame_descr *f = (frame_descr*)Field(callstacks[i], j);
      uint32 h = caml_hash_mix_uint32(caml_hash_mix_intnat(0, (intnat)f), 0);
      uintnat pos = j == 0 ? buffer[i].alloc_frame_pos : 0;
#else
      code_t c = (code_t)Field(callstacks[i], j);
      uint32 h = caml_hash_mix_intnat(0, (intnat)c);
#endif
      uintnat mask = (1LL << loc_htbl.lsize) - 1;
      uintnat k;

      for(k = h & mask; loc_htbl.htbl[k].isample != -1; k = (k+1) & mask) {
#ifdef NATIVE_CODE
        uintnat j_k = loc_htbl.htbl[k].jloc;
        uintnat isample = loc_htbl.htbl[k].isample;
        frame_descr *f_k = (frame_descr*)Field(callstacks[isample], j_k);
        uint32 pos_k = j_k == 0 ? buffer[isample].alloc_frame_pos : 0;

        if(f == f_k && pos == pos_k)
          break;
#else
        code_t c_k = (code_t)
          Field(callstacks[loc_htbl.htbl[k].isample], loc_htbl.htbl[k].jloc);
        if(c == c_k)
          break;
#endif
      }

      if(loc_htbl.htbl[k].isample == -1) {
#ifdef NATIVE_CODE
        caml_extract_location_info(f, pos, &li);
#else
        caml_extract_location_info(events, c, &li);
#endif
        Store_field(trace, j, value_of_locinfo(&li, h));

        loc_htbl.htbl[k].isample = i;
        loc_htbl.htbl[k].jloc = j;
        loc_htbl.htbl[k].h = h;
        loc_htbl.slots++;
        loc_htbl_resize(&loc_htbl);
        if(loc_htbl.htbl == NULL) goto oom;
      } else
        Store_field(trace, j,
          Field(Field(Field(res, loc_htbl.htbl[k].isample), 0), loc_htbl.htbl[k].jloc));
    }
  }

  End_roots();

  free(loc_htbl.htbl);
  free(buffer);
  free(callstacks);

  CAMLreturn(res);

  oom:
  free(buffer);
  free(callstacks);
  caml_raise_out_of_memory();
}
