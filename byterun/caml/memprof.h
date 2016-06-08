#ifndef CAML_MEMPROF_H
#define CAML_MEMPROF_H

#include "roots.h"

extern void caml_memprof_renew_minor_sample(void);

extern value* caml_memprof_young_limit;

extern value caml_memprof_track_alloc_shr(value block);
extern void caml_memprof_postpone_track_alloc_shr(value block);
extern void caml_memprof_track_young(uintnat wosize);
extern void caml_memprof_track_interned(header_t* block, header_t* blockend);
extern void caml_memprof_handle_postponed();

#ifdef NATIVE_CODE
extern void caml_memprof_call_gc_end(double exceeded_by);
extern double caml_memprof_call_gc_begin(void);
#endif

#endif /* CAML_MEMPROF_H */
