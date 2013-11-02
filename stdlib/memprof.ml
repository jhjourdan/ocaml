type ctrl = {
    lambda : float;
    dummy : unit
}

external get_ctrl : unit -> ctrl = "caml_memprof_get"
external set_ctrl : ctrl -> unit = "caml_memprof_set"

external clear : unit -> unit = "caml_memprof_clear"
external estimate_consumed_memory : unit -> int =
  "caml_memprof_estimate_consumed_memory"
