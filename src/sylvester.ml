open Owl
open Types

let print_dim x =
  let d1, d2 = Algodiff.D.Mat.shape x in
  Printf.printf "%i, %i\n" d1 d2

module Make (P : PT) = struct
  let dim = P.dim
  let flow_length = P.flow_length

  type layer =
    { mutable w : Algodiff.D.t
    ; mutable u : Algodiff.D.t
    ; mutable b : Algodiff.D.t
    ; mutable wp : Algodiff.D.t
    ; mutable up : Algodiff.D.t
    ; mutable bp : Algodiff.D.t }

  let prms =
    Array.init flow_length (fun _ ->
        { w = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim dim
        ; u = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim dim
        ; b = Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim 1
        ; wp = Algodiff.D.Mat.ones dim dim
        ; up = Algodiff.D.Mat.ones dim dim
        ; bp = Algodiff.D.Mat.ones dim 1 } )

  let reset_prms () =
    Array.iter
      (fun l ->
        l.w <- Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim dim;
        l.u <- Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim dim;
        l.b <- Algodiff.D.Mat.uniform ~a:(-0.01) ~b:0.01 dim 1;
        l.wp <- Algodiff.D.Mat.ones dim dim;
        l.up <- Algodiff.D.Mat.ones dim dim;
        l.bp <- Algodiff.D.Mat.ones dim 1 )
      prms

  let onestep z l =
    let open Algodiff.D.Maths in
    let w = l.w in
    let u = l.u in
    let b = l.b in
    let u = triu u in
    let q, r = Algodiff.D.Linalg.qr w in
    let a = (u *@ (transpose q) *@ z) + b in
    let x = tanh a in
    let z = z + (w *@ x) in
    let ld =
      let h = (F 1. - sqr x) in
      let d = transpose (diag u * diag r) in
      F 1. + (h * d) |> abs |> log
    in
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
