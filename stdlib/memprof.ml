type ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
}

external get_ctrl : unit -> ctrl = "caml_memprof_get"
external set_ctrl : ctrl -> unit = "caml_memprof_set"

external clear : unit -> unit = "caml_memprof_clear"

type loc_info = {
    li_hash : int;
    li_filename : string;
    li_line : int;
    li_start_chr : int;
    li_end_chr : int
}

type sample = {
    callstack : loc_info array;
    size : int;
    occurences : int;
}

external dump_samples : unit -> sample array = "caml_memprof_dump_samples"
