open Owl

let ring z =
  let open Algodiff.D.Maths in
  let z1 = get_row z 0 in
  let z2 = get_row z 1 in
  let nc =
    let nc = sqrt (sqr z1 + sqr z2) in
    (nc - F 4.) / F 0.4 |> sqr
  in
  let z1_ = (z1 - F 2.) / F 0.8 in
  let z2_ = (z1 + F 2.) / F 0.8 in
  let x1 = F 0.5 * sqr z1_ |> neg |> exp in
  let x2 = F 0.5 * sqr z2_ |> neg |> exp in
  (F 0.5 * nc) - log (x1 + x2)  

let ring2 z =
  let open Algodiff.D.Maths in
  let z1 = get_row z 0 in
  let z2 = get_row z 1 in
  let nc =
    let nc = sqrt (sqr z1 + sqr z2) in
    (nc - F 2.) / F 0.4 |> sqr
  in
  let z1_ = (z1 - F 2.) / F 0.6 in
  let z2_ = (z1 + F 2.) / F 0.6 in
  let x1 = F 0.5 * sqr z1_ |> neg |> exp in
  let x2 = F 0.5 * sqr z2_ |> neg |> exp in
  (F 0.5 * nc) - log (x1 + x2) 

let wave z =
  let open Algodiff.D.Maths in
  let z1 = get_row z 0 in
  let z2 = get_row z 1 in
  let w = (F (Const.pi *. 0.5)) * z1 |> sin in
  let x = (z2 - w) / (F 0.4) |> sqr in
  (F 0.5) * x 




