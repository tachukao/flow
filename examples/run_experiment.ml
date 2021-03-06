open Owl
open Flow

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let in_dir = Printf.sprintf "%s/%s" dir

let () =
  let n_samples = 5000 in
  let open Algodiff.D in
  let module F =
    Sylvester.Make (struct
      let flow_length = 5
      let dim = 2
    end)
  in
  F.reset_prms ();
  train
    ~eta:0.001
    ~max_iter:10000
    ~print_every:1000
    ~batch_size:40
    (module F)
    Potentials.ring;
  let samples = sample (module F) n_samples in
  Owl.Mat.(save_txt (transpose (unpack_arr samples)) ~out:(in_dir "samples"))
