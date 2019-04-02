open Owl

let () =
  let open Algodiff.D in
  let module F = Sylvester.Make (struct
    let flow_length = 32
  end) in
  Flow.train (module F) Density.wave;
  let z = Mat.gaussian 2 5000 in
  let samples, _ = Flow.flow (module F) z in
  Owl.Mat.(save_txt (transpose (unpack_arr samples)) "results/samples")
