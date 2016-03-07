type ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
}

val get_ctrl : unit -> ctrl
val set_ctrl : ctrl -> unit

val reset : unit -> unit

type sample = {
    size : int;
    occurences : int;
    callstack : Printexc.raw_backtrace
}
