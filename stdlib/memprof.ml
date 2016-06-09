type callback_kind =
  | Minor
  | Major
  | Major_postponed
  | Serialized

type 'a callback = callback_kind -> int -> int -> Printexc.raw_backtrace ->
                   (Obj.t, 'a) Ephemeron.K1.t option

type 'a ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
    callback : 'a callback
}

let stopped_ctrl = {
    lambda = 0.; dumpped_callstack_size = 0;
    callback = fun _ _ _ _ -> assert false
}

external set_ctrl : 'a ctrl -> unit = "caml_memprof_set"

let start = set_ctrl
let stop () = set_ctrl stopped_ctrl
