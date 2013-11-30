/***********************************************************************/
/*                                                                     */
/*                                OCaml                                */
/*                                                                     */
/*            Xavier Leroy, projet Cristal, INRIA Rocquencourt         */
/*                                                                     */
/*  Copyright 2001 Institut National de Recherche en Informatique et   */
/*  en Automatique.  All rights reserved.  This file is distributed    */
/*  under the terms of the GNU Library General Public License, with    */
/*  the special exception on linking described in file ../LICENSE.     */
/*                                                                     */
/***********************************************************************/

#ifndef CAML_BACKTRACE_H
#define CAML_BACKTRACE_H

#include "mlvalues.h"
#ifdef NATIVE_CODE
#include "../asmrun/stack.h"
#endif

CAMLextern int caml_backtrace_active;
CAMLextern int caml_backtrace_pos;
CAMLextern code_t * caml_backtrace_buffer;
CAMLextern value caml_backtrace_last_exn;
CAMLextern char * caml_cds_file;

CAMLprim value caml_record_backtrace(value vflag);
#ifndef NATIVE_CODE
extern void caml_stash_backtrace(value exn, code_t pc, value * sp, int reraise);
#endif
CAMLextern void caml_print_exception_backtrace(void);

#ifdef NATIVE_CODE
extern frame_descr * caml_next_frame_descriptor(uintnat * pc, char ** sp);
#else
extern code_t caml_next_frame_pointer(value ** sp, value ** trapsp);
#endif

value caml_get_current_callstack(value max_frames_value);

struct loc_info {
  int loc_valid;
  int loc_is_raise;
  char * loc_filename;
  int loc_lnum;
  int loc_startchr;
  int loc_endchr;
};

#ifdef NATIVE_CODE
/* Extract location information for the given frame descriptor */
/* alloc_id is the identifier of the allocation in a comballoc group */

void caml_extract_location_info(frame_descr * d, intnat alloc_id,
                                /*out*/ struct loc_info * li);

#else
/* Extract location information for the given PC */

void caml_extract_location_info(value events, code_t pc,
                                /*out*/ struct loc_info * li);

CAMLextern char *caml_read_debug_info_error;
value caml_read_debug_info(void);
#endif

#ifdef NATIVE_CODE
typedef struct {
  uintnat whsize;
  uint32 loc1, loc2;
} alloc_descr;
#endif

#endif /* CAML_BACKTRACE_H */
