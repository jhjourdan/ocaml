#include <math.h>
#include <string.h>
#include "caml/memprof.h"
#include "caml/backtrace_prim.h"
#include "caml/fail.h"
#include "caml/alloc.h"
#include "caml/callback.h"
#include "caml/weak.h"
#include "caml/signals.h"
#ifdef NATIVE_CODE
#include "stack.h"
#endif

static uint32_t mt_state[624];
static uint32_t mt_index;

/* [lambda] is the mean number of samples for each allocated word (including
   block headers). */
static double lambda = 0;
static double lambda_rec = INFINITY;
static intnat callstack_size = 0;
static value memprof_callback = Val_unit;

/* Position of the next sample in the minor heap. Equals
   [caml_young_alloc_start] if no sampling is planned in the current
   minor heap. */
value* caml_memprof_young_limit;

/* The continuous position of the next minor sample within the sampled
   word. Always in [0, 1[. */
static double next_sample_round;

/* Whether memprof has been initialized.  */
static int init = 0;

/**** Statistical sampling ****/

/* Taken from :
   https://github.com/jhjourdan/SIMD-math-prims
*/
inline static float logapprox(float val) {
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

inline static float expapprox(float val) {
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

/* Max returned value : 2^30-2. Assumes lambda >= 0. */
static int32_t mt_generate_poisson(double lambda) {
  Assert(lambda >= 0 && !isinf(lambda));

  if(lambda == 0)
    return 0;

  if(lambda < 20) {
    double p;
    int32_t k;
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

/**** Interface with the OCaml code. ****/

void set_lambda(double l) {
  lambda = l;
  lambda_rec = l == 0 ? INFINITY : 1/l;
  caml_memprof_renew_minor_sample();
}

CAMLprim value caml_memprof_set(value v) {
  CAMLparam1(v);
  double l = Double_val(Field(v, 0));
  intnat sz = Long_val(Field(v, 1));

  if(sz < 0 || !(l >= 0.) || l > 1.)
    caml_failwith("caml_memprof_set");

  set_lambda(l);
  callstack_size = sz;
  memprof_callback = Field(v, 2);

  if(!init) {
    int i;
    mt_index = 624;
    mt_state[0] = 42;
    for(i = 1; i < 624; i++)
      mt_state[i] = 0x6c078965 * (mt_state[i-1] ^ (mt_state[i-1] >> 30)) + i;

    caml_register_global_root(&memprof_callback);
  }

  CAMLreturn(Val_unit);
}

/* Cf. Memprof.callstack_kind */
enum ml_callback_kind {
  Minor = Val_long(0),
  Major = Val_long(1),
  Major_postponed = Val_long(2),
  Serialized = Val_long(3)
};

static value do_callback(intnat wosize, int32_t occurences, value callstack,
                         enum ml_callback_kind cb_kind) {
  Assert(occurences > 0);
  value args[4] =
    { cb_kind, Val_long(wosize), Val_long(occurences), callstack };
  return caml_callbackN_exn(memprof_callback, 4, args);
}

/**** Sampling procedures ****/

/* Shifts the next sample in the minor heap by [n] words. Essentially,
   this tells the sampler to ignore the next [n] words of the minor
   heap. */
static void shift_sample(uintnat n) {
  if(caml_memprof_young_limit - caml_young_alloc_start > n)
    caml_memprof_young_limit -= n;
  else
    caml_memprof_young_limit = caml_young_alloc_start;
  caml_update_young_limit();
}

/* Renew the next sample in the minor heap. This needs to be called
   after each minor sampling and after each minor collection. In
   practice, because we disable sampling during callbacks, this is
   called at each sampling (including major ones). These extra calls
   do not change the statistical properties of the sampling because of
   the memorylessness of the exponential distribution. */
void caml_memprof_renew_minor_sample(void) {
  double exp = mt_generate_exponential();
  uintnat max = caml_young_ptr - caml_young_alloc_start, exp_int;
  if(exp < max  && (exp_int = (uintnat)exp) < max) {
    caml_memprof_young_limit = caml_young_ptr - exp_int;
    next_sample_round = exp - exp_int;
    Assert(0 <= next_sample_round && next_sample_round < 1);
  } else
    caml_memprof_young_limit = caml_young_alloc_start;
  caml_update_young_limit();
}

static value capture_callstack(int avoid_gc) {
  return caml_get_current_callstack_impl(callstack_size, avoid_gc);
}

/* Called when exceeding the threshold for the next sample in the
   minor heap, from the C code (the handling is different when call
   from natively compiled OCaml code). */
void caml_memprof_track_young(uintnat wosize) {
  CAMLparam0();
  CAMLlocal2(ephe, callstack);
  uintnat whsize = Whsize_wosize(wosize);
  double rest =
    (caml_memprof_young_limit - caml_young_ptr) - next_sample_round;

  Assert(rest > 0 && lambda > 0);

  int32_t occurences = mt_generate_poisson(rest*lambda) + 1;

  caml_young_ptr += whsize;
  //We should not allocate before this point

  caml_memprof_handle_postponed();

  double old_lambda = lambda;
  set_lambda(0);
  callstack = capture_callstack(0);
  ephe = do_callback(wosize, occurences, callstack, Minor);
  set_lambda(old_lambda); // Calls [caml_memprof_renew_minor_sample]
  if (Is_exception_result(ephe)) caml_raise(Extract_exception(ephe));

  if(caml_young_ptr - whsize < caml_young_trigger)
    caml_gc_dispatch();

  // We should not allocate after this point
  caml_young_ptr -= whsize;
  shift_sample(whsize); // Make sure this block is not going th be
                        // sampled again.
  if(Is_block(ephe))
    caml_ephe_set_key(Field(ephe, 0), Val_long(0), Val_hp(caml_young_ptr));

  CAMLreturn0;
}

/* Called when allocating in the major heap, except when allocating
   during a minor collection or for a serialization, or when we do not
   want to call the GC. */
value caml_memprof_track_alloc_shr(value block) {
  CAMLparam1(block);
  CAMLlocal2(ephe, callstack);

  Assert(Is_in_heap(block));

  int32_t occurences = mt_generate_poisson(lambda * Whsize_val(block));

  caml_memprof_handle_postponed();

  if(occurences > 0) {
    double old_lambda = lambda;
    set_lambda(0);
    callstack = capture_callstack(0);
    ephe = do_callback(Wosize_val(block), occurences, callstack, Major);
    set_lambda(old_lambda);
    if (Is_exception_result(ephe)) caml_raise(Extract_exception(ephe));

    if(Is_block(ephe))
      caml_ephe_set_key(Field(ephe, 0), Val_long(0), block);
  }

  CAMLreturn (block);
}

struct postponed_block {
  value block;
  value callstack;
  int32_t occurences;
  struct postponed_block* next;
} *postponed_head = NULL;

/* If [caml_alloc_shr] is called in a context where calling the GC
   should be avoided (historically, alloc_shr does not call the GC),
   we postpone the call of the callback. This function is analoguous
   to [caml_memprof_track_alloc_shr], except that if anything is
   sampled, it stores the block in the todo-list so that the callback
   call is performed later. */
void caml_memprof_postpone_track_alloc_shr(value block) {
  int32_t occurences = mt_generate_poisson(lambda * Whsize_val(block));
  Assert(Is_in_heap(block));
  if(occurences > 0) {
    struct postponed_block* pb =
      (struct postponed_block*)malloc(sizeof(struct postponed_block));
    if(pb == NULL) return;
    pb->block = block;
    caml_register_generational_global_root(&pb->block);
    double old_lambda = lambda;
    set_lambda(0);
    pb->callstack = capture_callstack(1);
    set_lambda(old_lambda);
    caml_register_generational_global_root(&pb->callstack);
    pb->occurences = occurences;
    pb->next = postponed_head;
    postponed_head = pb;
#ifndef NATIVE_CODE
    caml_something_to_do = 1;
#else
    caml_young_limit = caml_young_alloc_end;
#endif
  }
}

void caml_memprof_handle_postponed() {
  struct postponed_block *p, *q;
  double old_lambda = lambda;
  value ephe;

  if(postponed_head == NULL)
    return;

  // We first revert the list
  p = postponed_head;
  q = postponed_head->next;
  p->next = NULL;
  while(q != NULL) {
    struct postponed_block* next = q->next;
    q->next = p;
    p = q;
    q = next;
  }
  postponed_head = NULL;

#define NEXT_P \
  { caml_remove_generational_global_root(&p->callstack); \
    caml_remove_generational_global_root(&p->block);     \
    struct postponed_block* next = p->next;              \
    free(p);                                             \
    p = next; }

  set_lambda(0);
  // We then do the actual iteration on postponed blocks
  while(p != NULL) {
    ephe = do_callback(Wosize_val(p->block), p->occurences, p->callstack,
                       Major_postponed);
    if (Is_exception_result(ephe)) {
      set_lambda(old_lambda);
      // In the case of an exception, we just forget the entire list.
      while(p != NULL) NEXT_P;
      caml_raise(Extract_exception(ephe));
    }
    if(Is_block(ephe))
      caml_ephe_set_key(Field(ephe, 0), Val_long(0), p->block);
    NEXT_P;
  }
  set_lambda(old_lambda);
}

void caml_memprof_track_interned(header_t* block, header_t* blockend) {
  CAMLparam0();
  CAMLlocal2(ephe, callstack);

  value* sampled = NULL;
  int32_t* occurences = NULL;
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
      occurences = (int32_t*)malloc(sizeof(int32_t) * sz);
      if(occurences == NULL) {
        free(sampled);
        CAMLreturn0;
      }
    } else if(j >= sz) {
      value* newsampled;
      int32_t* newoccurences;
      sz *= 2;
      newsampled = (value*)realloc(sampled, sizeof(value) * sz);
      if(newsampled == NULL) {
        free(sampled); free(occurences);
        CAMLreturn0;
      }
      sampled = newsampled;
      newoccurences = (int32_t*)realloc(occurences, sizeof(int32_t) * sz);
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

  // We should not allocate before this point
  CAMLxparamN(sampled, j);

  caml_memprof_handle_postponed();

  double old_lambda = lambda;
  set_lambda(0);
  callstack = capture_callstack(0);
  for(i = 0; i < j; i++) {
    ephe = do_callback(Wosize_val(sampled[i]), occurences[i], callstack,
                       Serialized);
    if (Is_exception_result(ephe)) {
      free(sampled);
      free(occurences);
      set_lambda(old_lambda);
      caml_raise(Extract_exception(ephe));
    }
    if(Is_block(ephe))
      caml_ephe_set_key(Field(ephe, 0), Val_long(0), sampled[i]);
  }

  free(sampled);
  free(occurences);
  set_lambda(old_lambda);
  CAMLreturn0;
}

#ifdef NATIVE_CODE

/* Sampling in the minor heap for native code is different from
   sampling for C-allocated blocks, because the compiler back-end
   merge allocations for better performance. Therefore, we need to
   revert this mechanism and recover the pointers to the different
   blocks being alocated. */

/* Allocation in compiled OCaml code proceeds by retrying until
   [caml_young_ptr] is under [caml_young_limit]. If there are several
   tries, we need to make sure that the block is sampled only once. To
   that end:
     - The sampling is ignored if the allocation fails because of the
       minor heap being full (or if a gc minor/major collection has
       been requested).

     - As soon as the sampling is proceeded (i.e., we enter
       [caml_memprof_call_gc_end]), we make sure there is enough room
       in the minor heap and we shift the threshold so that the block
       is not going to be sampled. Note that because a signal can (in
       theory) be raised just before the allocation, this is not
       guaranteed that the allocation will actually take place even
       with these precautions, but we ignore this eventuallity.

   See [caml_garbage_collection] in signals_asm.c.
 */

/* If the allocation has failed because the block has to be sampled,
   we record by how much the threshold has been exceeded. */
double caml_memprof_call_gc_begin(void) {
  intnat exceeded_int = caml_memprof_young_limit - caml_young_ptr;
  if(exceeded_int > 0) {
    double exceeded = exceeded_int - next_sample_round;
    caml_memprof_renew_minor_sample();
    return exceeded;
  } else
    return 0;
}

void caml_memprof_call_gc_end(double next_sample) {
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

  if(next_sample > 0) {
    Assert(next_sample <= tot_whsize);

    uintnat offs = 0;
    struct smp { uintnat offs, sz, alloc_frame_pos; int32_t occurences; }
      *samples, *p;
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
          p = (struct smp*)realloc(samples, sz*sizeof(struct smp));
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

    caml_memprof_handle_postponed();

    CAMLlocalN(ephes, n_samples);
    double old_lambda = lambda;
    set_lambda(0);
    callstack = capture_callstack(0);
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

      ephes[i] =
        do_callback(samples[i].sz, samples[i].occurences, callstack_cur, Minor);
      if(Is_exception_result(ephes[i])) {
        free(samples);
        set_lambda(old_lambda);
        caml_raise(Extract_exception(ephes[i]));
      }
    }
    set_lambda(old_lambda); // Calls caml_memprof_renew_minor_sample

    // Do everything we can to make sure that the allocation will take
    // place.
    if(caml_requested_major_slice || caml_requested_minor_gc ||
       caml_young_ptr - caml_young_trigger < tot_whsize)
      caml_gc_dispatch();
    // We should not allocate after this point

    for(i = 0; i < n_samples; i++) {
      value v = Val_hp(caml_young_ptr - samples[i].offs);

      /* It is important not to allocate anything or call the GC between
       * this point and the re-execution of the allocation code in assembly,
       * because blk.block would then be incorrect. */
      if(Is_block(ephes[i]))
        caml_ephe_set_key(Field(ephes[i], 0), Val_long(0), v);

      /* It is *not* garanteed that this block will actually be used,
       * because a signal can interrupt the allocation in the assembly
       * code. Thus, we put an abstract header on this unitialized
       * block. We do not consider this eventuality in statistical
       * computations. */
      Hd_val(v) = Make_header(samples[i].sz, Abstract_tag, Caml_black);
    }

    free(samples);
  } else if(caml_requested_major_slice || caml_requested_minor_gc ||
            caml_young_ptr - tot_whsize < caml_young_trigger)
    caml_gc_dispatch();

 abort:
  /* We prevent the next allocation to be sampled, as it already had its
   * chance in this call. */
  shift_sample(tot_whsize);

  CAMLreturn0;
}
#endif
