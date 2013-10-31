#ifndef CAML_MEMPROF_H
#define CAML_MEMPROF_H

#include "roots.h"

struct caml_memprof_tracked_block {
  value block;
#ifdef NATIVE_CODE
  uint32 loc1, loc2;
#else
  code_t loc;
#endif
  uint32 occurences;
};

extern struct caml_memprof_tracked_block* caml_memprof_tracked_blocks;
extern uintnat caml_memprof_tracked_blocks_end;

extern void caml_memprof_minor_gc(void);
extern void caml_memprof_major_gc(void);
extern void caml_memprof_do_weak_roots(scanning_action f);
extern void caml_memprof_reinit(void);

extern char* caml_memprof_young_limit;

extern void caml_memprof_track_one(value block, uintnat wosize);

#ifdef NATIVE_CODE
extern double caml_memprof_call_gc_begin_hook(void);
extern void caml_memprof_call_gc_end_hook(double);
#endif

extern double caml_memprof_get_lambda(void);
extern void caml_memprof_set_lambda(double);

#endif /* CAML_MEMPROF_H */
