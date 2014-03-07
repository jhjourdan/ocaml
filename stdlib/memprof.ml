type ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
}

external set_lambda : float -> unit = "caml_memprof_set"

let ctrl = ref { lambda = 0.; dumpped_callstack_size = 3 }
let get_ctrl () = !ctrl
let set_ctrl c =
  if c.dumpped_callstack_size < 0 ||
     not (c.lambda >= 0.) || c.lambda > 1. then
    raise (Invalid_argument "caml_memprof_set");
  set_lambda c.lambda;
  ctrl := c

external reset : unit -> unit = "caml_memprof_clear"

type sample = {
    callstack : Printexc.raw_backtrace;
    size : int;
    occurences : int;
}

external dump_samples : unit -> sample array = "caml_memprof_dump_samples"

external register_callback : (unit -> Printexc.raw_backtrace) -> unit
    = "caml_memprof_register_callback"
let _ = register_callback
    (fun () -> Printexc.get_callstack (!ctrl.dumpped_callstack_size+1))
