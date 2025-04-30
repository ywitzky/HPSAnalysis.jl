# Requirements for HPSAnalysis.jl

Julia installation after cloning from github:
```julia
    ] dev /path/to/HPSAnalysis.jl  
```

## (Test of) Setup Routine

The module named Setup allows to set up Calvados2 and Calvados3 via the Python wrapper of the simulation package named HOOMD. Therefore, one needs to manually install a python environment with the following packages

* Polyply
* martinize2 (Calvados3)
* HOOMD
* GSD

For the tests to work, one needs to activate the Python environment before starting Julia and safe the Path to the Python executable in the file "/HPSAnalysis/data/EnvironmentPath.txt" . 