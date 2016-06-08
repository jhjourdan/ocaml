open Ephemeron

type ctrl = {
    lambda : float;
    dumpped_callstack_size : int;
}

external set_ctrl : ctrl -> unit = "caml_memprof_set"
external get_ctrl : unit -> ctrl = "caml_memprof_get"

external register_callback :
  (int -> int -> Printexc.raw_backtrace -> (Obj.t, 'a) K1.t) -> unit
  = "caml_memprof_register_callback"

type sample = {
    size : int;
    occurences : int;
    callstack : Printexc.raw_backtrace
}

let min_buf_size = 1024
let empty_ephe = K1.create ()
let samples = ref (Array.make min_buf_size empty_ephe)
let n_samples = ref 0

let reset () =
  n_samples := 0;
  samples := [| empty_ephe |]

let clean () =
  let s = !samples in
  let sz = !n_samples in
  let rec aux i j =
    if i >= sz then j
    else if K1.check_key s.(i) then (s.(j) <- s.(i); aux (i+1) (j+1))
    else aux (i+1) j
  in
  let new_sz = aux 0 0 in
  Array.fill s new_sz (sz - new_sz) empty_ephe;
  n_samples := new_sz;
  if 8 * !n_samples <= Array.length s &&
     Array.length s > min_buf_size then
    samples := Array.sub s 0 (max min_buf_size (2 * !n_samples))

let push e =
  while !n_samples = Array.length !samples do
    clean ();
    if !n_samples = Array.length !samples then begin
      let rec alloc () =
        let res = Array.make (!n_samples * 2) empty_ephe in
        if Array.length res <> !n_samples * 2 then alloc ()
        else res
      in
      let s = alloc () in
      Array.blit !samples 0 s 0 !n_samples;
      samples := s
    end
  done;
  !samples.(!n_samples) <- e;
  incr n_samples

let _ =
  register_callback
    (fun size occurences callstack ->
       let ephe = K1.create () in
       K1.set_data ephe { callstack; size; occurences };
       push ephe;
       ephe)
