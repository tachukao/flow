open Owl

module type PT = sig
  val flow_length : int
end

module type FlowT = sig
  type layer

  val flow_length : int
  val prms : layer array
  val onestep : Algodiff.D.t -> layer -> Algodiff.D.t * Algodiff.D.t
  val tag_layer : int -> layer -> unit
  val update_layer : float -> layer -> unit
  val reset_prms : unit -> unit
end
