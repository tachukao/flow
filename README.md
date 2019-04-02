# Implementation of normalizing flows

### Usage
```ocaml
let samples = 
   let module F = Planar.Make (struct let flow_length = flow_length end) in
   Flow.train ~eta ~max_iter ~batch_size (module F) Potentials.ring;
   Flow.sample (module F) n_samples
```

### Example
```sh
mkdir results
dune exec examples/run_experiment.exe
```
