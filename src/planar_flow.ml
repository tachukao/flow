open Owl

let safe_log z = Algodiff.D.(Maths.(log (z + F 1E-9)))
let batch_size = 40
let flow_length = 32
let max_iter = 20000
let eta = Algodiff.D.pack_flt 0.001
let density = Density.wave

type layer =
  { mutable w : Algodiff.D.t
  ; mutable b : Algodiff.D.t
  ; mutable u : Algodiff.D.t
  ; mutable wp : Algodiff.D.t
  ; mutable bp : Algodiff.D.t
  ; mutable up : Algodiff.D.t }

let prms =
  Array.init flow_length (fun _ ->
      { w = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 2 1
      ; b = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 1 1
      ; u = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 2 1
      ; wp = Algodiff.D.Mat.ones 2 1
      ; bp = Algodiff.D.Mat.ones 1 1
      ; up = Algodiff.D.Mat.ones 2 1 } )


let onestep z w b u =
  let open Algodiff.D.Maths in
  let a = (transpose w *@ z) + b in
  let r = tanh a in
  let z = z + (r * u) in
  let psi = (F 1. - sqr r) * w in
  let ld = F 1. + (transpose u *@ psi) |> abs |> safe_log in
  z, ld

let rec planar_flow i z lj =
  if i < flow_length
  then
    let w = prms.(i).w in
    let u = prms.(i).u in
    let b = prms.(i).b in
    let z, ld = onestep z w b u in
    let lj = Algodiff.D.Maths.(lj + ld) in
    planar_flow (succ i) z lj
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
      let z, lj = planar_flow 0 z (F 0.) in
      Maths.(sum' (neg (safe_log (density z)) - lj))
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
  let samples, _ = planar_flow 0 z (F 0.) in
  Owl.Mat.(save_txt (transpose (unpack_arr samples)) "results/samples")



