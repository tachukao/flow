open Owl

let batch_size = 40
let flow_length = 16
let max_iter = 20000
let eta = Algodiff.D.pack_flt 0.001
let density = Density.wave

type layer =
  { mutable w : Algodiff.D.t
  ; mutable u : Algodiff.D.t
  ; mutable b : Algodiff.D.t
  ; mutable wp : Algodiff.D.t
  ; mutable up : Algodiff.D.t
  ; mutable bp : Algodiff.D.t }

let prms =
  Array.init flow_length (fun _ ->
      { w = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 2 2
      ; u = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 2 2
      ; b = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 2 1
      ; wp = Algodiff.D.Mat.ones 2 2
      ; up = Algodiff.D.Mat.ones 2 2
      ; bp = Algodiff.D.Mat.ones 2 1 } )

let onestep z l =
  let open Algodiff.D.Maths in
  let w = l.w in
  let u = l.w in
  let b = l.w in
  let q, r = qr w in
  let u = triu u in
  let a = (u *@ transpose q *@ z) + b in
  let x = tanh a in
  let z = z + (w *@ x) in
  let ld =
    let h = F 1. - sqr x in
    let d = transpose (diag u * diag r) in
    (h * d) + F 1. |> abs |> log
  in
  z, ld

let rec flow i z lj =
  if i < flow_length
  then
    let l = prms.(i) in
    let z, ld = onestep z l in
    let lj = Algodiff.D.Maths.(lj + ld) in
    flow (succ i) z lj
  else z, lj

let rec optimize iter loss =
  let open Algodiff.D in
  if iter < max_iter && loss > F (-1.)
  then (
    (*let eta = F (eta *. Owl.Maths.(exp (-0.999 *. float iter))) in*)
    let z = Mat.gaussian 2 batch_size in
    let t = tag () in
    Array.iter
      (fun l ->
        l.w <- make_reverse l.w t;
        l.u <- make_reverse l.u t;
        l.b <- make_reverse l.b t )
      prms;
    let loss =
      let z, lj = flow 0 z (F 0.) in
      Maths.(sum' (neg (log (density z)) - lj))
    in
    reverse_prop (F 1.) loss;
    Array.iter
      (fun l ->
        let dw = adjval l.w |> primal in
        let du = adjval l.u |> primal in
        let db = adjval l.b |> primal in
        l.w <- Maths.(primal l.w - (eta * dw / sqrt l.wp)) |> primal;
        l.u <- Maths.(primal l.u - (eta * du / sqrt l.up)) |> primal;
        l.b <- Maths.(primal l.b - (eta * db / sqrt l.bp)) |> primal;
        l.wp <- Maths.((F 0.9 * l.wp) + (F 0.1 * sqr dw)) |> primal;
        l.up <- Maths.((F 0.9 * l.up) + (F 0.1 * sqr du)) |> primal;
        l.bp <- Maths.((F 0.9 * l.bp) + (F 0.1 * sqr db)) |> primal )
      prms;
    if iter mod 1000 = 0
    then Printf.printf "iter: %i | cost %f | lr %f %!\n" iter (unpack_flt loss) (unpack_flt eta);
    optimize (succ iter) loss )

let () =
  let open Algodiff.D in
  optimize 0 (F 1E9);
  let z = Mat.gaussian 2 5000 in
  let samples, _ = flow 0 z (F 0.) in
  Owl.Mat.(save_txt (transpose (unpack_arr samples)) "results/samples")
