type ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
}

external get_ctrl : unit -> ctrl = "caml_memprof_get"
external set_ctrl : ctrl -> unit = "caml_memprof_set"

external clear : unit -> unit = "caml_memprof_clear"

type loc_info =
    string  (* filename *)
      * int (* line number *)
      * int (* start char *)
      * int (* end char *)

type sample = {
    callstack : loc_info array;
    size : int;
    occurences : int;
}

external dump_samples : unit -> sample array = "caml_memprof_dump_samples"
