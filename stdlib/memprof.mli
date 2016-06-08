type callback_kind =
  | Minor
  | Major
  | Major_postponed
  | Serialized

type 'a ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
    callback :
      callback_kind -> int -> int -> Printexc.raw_backtrace ->
      (Obj.t, 'a) Ephemeron.K1.t
}

val start : 'a ctrl -> unit
val stop : unit -> unit
