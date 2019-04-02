open Owl


let () =
  let open Algodiff.D in
  let module F = Sylvester.Make (struct
    let flow_length = 5
  end) in
  F.reset_prms ();
  Flow.train ~eta:0.001 ~max_iter:10000 ~print_every:1000 ~batch_size:40 (module F) Potentials.ring;
  let z = Mat.gaussian 2 5000 in
  let samples, _ = Flow.flow (module F) z in
  Owl.Mat.(save_txt (transpose (unpack_arr samples)) "results/samples")
