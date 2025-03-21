Source code used for the master thesis "Numerical methods for localizing flow in the geosciences", handed in on March 21st 2025.

In the following a short explanation of the repository structure:

- `bTimes/` contains log files from the benchmarks that are used in some figure plots. 
- `Figures/` contains all the scripts for plotting the figures from the thesis as well as the raw .svg files for the illustrations. In all the `.jl` scripts there is a bool `p2f` which if enabled saves a `.png` file to the folder.  
- `Benchmark.jl` defines some routines to benchmark the solver with various options (e.g. differnt grid spacings). If these functions are executed they will overwrite the data in `bTimes/`.
- `benchRapidLoc.jl` contains the solver with a parameter set that results und rapid shear localization and is used by the routines defined in `Benchmark.jl`
- `ShearLocalization.jl` is the raw main solver with an arbitrary paramter set which can be modified.
- `helperFunctions.jl` contains some smaller routines e.g. a function to compute the minimum timestep.   

Created with Julia version 1.10.8
