open Owl
open Types

module Make (P : PT) = struct
  let dim = P.dim
  let flow_length = P.flow_length

  type layer =
    { mutable w : Algodiff.D.t
    ; mutable b : Algodiff.D.t
    ; mutable u : Algodiff.D.t
    ; mutable wp : Algodiff.D.t
    ; mutable bp : Algodiff.D.t
    ; mutable up : Algodiff.D.t }

  let prms =
    Array.init flow_length (fun _ ->
        { w = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim 1
        ; b = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 1 1
        ; u = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim 1
        ; wp = Algodiff.D.Mat.ones dim 1
        ; bp = Algodiff.D.Mat.ones 1 1
        ; up = Algodiff.D.Mat.ones dim 1 } )

  let reset_prms () =
    Array.iter
      (fun l ->
        l.w <- Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim 1;
        l.u <- Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim 1;
        l.b <- Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 1 1;
        l.wp <- Algodiff.D.Mat.ones dim 1;
        l.up <- Algodiff.D.Mat.ones dim 1;
        l.bp <- Algodiff.D.Mat.ones 1 1 )
      prms

  let onestep z l =
    let open Algodiff.D.Maths in
    let w = l.w in
    let b = l.b in
    let u = l.u in
    let a = (transpose w *@ z) + b in
    let r = tanh a in
    let z = z + (r * u) in
    let psi = (F 1. - sqr r) * w in
    let ld = F 1. + (transpose u *@ psi) |> abs |> log in
    z, ld

  let tag_layer t l =
    l.w <- Algodiff.D.make_reverse l.w t;
    l.u <- Algodiff.D.make_reverse l.u t;
    l.b <- Algodiff.D.make_reverse l.b t

  let update_layer eta =
    let open Algodiff.D in
    let eta = F eta in
    fun l ->
      let dw = adjval l.w |> primal in
      let du = adjval l.u |> primal in
      let db = adjval l.b |> primal in
      l.w <- Maths.(primal l.w - (eta * dw / sqrt l.wp)) |> primal;
      l.u <- Maths.(primal l.u - (eta * du / sqrt l.up)) |> primal;
      l.b <- Maths.(primal l.b - (eta * db / sqrt l.bp)) |> primal;
      l.wp <- Maths.((F 0.9 * l.wp) + (F 0.1 * sqr dw)) |> primal;
      l.up <- Maths.((F 0.9 * l.up) + (F 0.1 * sqr du)) |> primal;
      l.bp <- Maths.((F 0.9 * l.bp) + (F 0.1 * sqr db)) |> primal
end
