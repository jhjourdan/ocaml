type ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
}

val get_ctrl : unit -> ctrl
val set_ctrl : ctrl -> unit

val reset : unit -> unit

type sample = {
    callstack : Printexc.raw_backtrace;
    size : int;
    occurences : int;
}

val dump_samples : unit -> sample array
