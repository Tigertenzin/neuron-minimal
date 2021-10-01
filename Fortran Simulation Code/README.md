# Fortran Simulation Code 

Different algorithms for simulating the minimal neuron model. 

* The `ER` version directly simulates the model as described (step by step). 
* The `randWalk` version maps the model's dynamics onto a random walker problem (where a given avalanche is a single random walker, and there is limited space that can be occupied by these random walkers). 
* The `histograms.f90` takes in the raw simulation data and returns log-binned histograms. 