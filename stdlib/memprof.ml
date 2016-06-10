type alloc_kind =
  | Minor
  | Major
  | Major_postponed
  | Serialized

type sample_info = {
    n_samples: int; kind: alloc_kind; tag: int;
    size: int; callstack: Printexc.raw_backtrace;
}

type 'a callback = sample_info -> (Obj.t, 'a) Ephemeron.K1.t option

type 'a ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
    callback : 'a callback
}

let stopped_ctrl = {
    lambda = 0.; dumpped_callstack_size = 0;
    callback = fun _ -> assert false
}

external set_ctrl : 'a ctrl -> unit = "caml_memprof_set"

let start = set_ctrl
let stop () = set_ctrl stopped_ctrl
