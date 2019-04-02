open Owl
open Types

let flow (module F : FlowT) z =
  let rec run i z lj =
    if i < F.flow_length
    then
      let l = F.prms.(i) in
      let z, ld = F.onestep z l in
      let lj = Algodiff.D.Maths.(lj + ld) in
      run (succ i) z lj
    else z, lj
  in
  run 0 z (F 0.)

let train ?(max_iter = 1000) ?(batch_size = 10) ?(eta = 1.) (module F : FlowT) density =
  let rec optimize iter loss =
    let open Algodiff.D in
    if iter < max_iter && loss > F (-1.)
    then (
      let z = Mat.gaussian 2 batch_size in
      let t = tag () in
      Array.iter (fun l -> F.tag_layer t l) F.prms;
      let loss =
        let z, lj = flow (module F) z in
        Maths.(sum' (neg (log (density z)) - lj))
      in
      reverse_prop (F 1.) loss;
      Array.iter (F.update_layer eta) F.prms;
      if iter mod 1000 = 0 then Printf.printf "iter: %i | cost %f%!\n" iter (unpack_flt loss);
      optimize (succ iter) loss )
  in
  optimize 0 (F 1E9)
