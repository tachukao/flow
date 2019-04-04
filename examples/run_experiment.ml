open Owl


let () =
  let n_samples = 5000 in
  let open Algodiff.D in
  let module F = Sylvester.Make (struct
    let flow_length = 5
    let dim = 2
  end) in
  F.reset_prms ();
  Flow.train ~eta:0.001 ~max_iter:10000 ~print_every:1000 ~batch_size:40 (module F) Potentials.ring;
  let samples = Flow.sample (module F) n_samples in
  Owl.Mat.(save_txt (transpose (unpack_arr samples)) "results/samples")
