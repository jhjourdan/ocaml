type ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
}

val get_ctrl : unit -> ctrl
val set_ctrl : ctrl -> unit

val clear : unit -> unit

type loc_info = {
    li_hash : int;
    li_filename : string;
    li_line : int;
    li_start_chr : int;
    li_end_chr : int;
}
type sample = {
    callstack : loc_info array;
    size : int;
    occurences : int;
}

val dump_samples : unit -> sample array
