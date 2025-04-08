# InhomogeneousBathymetryEuler

This code computes the incompressible Euler equations for an irrotational fluid in the form of potential flow with a linearized surface condition with bathymetry. It should be executable by calling 
```julia
julia> using Pkg

julia> Pkg.instantiate()
```
in the main directory. 
The structure of the file system is as follows:
* ``src/`` contains the source code
* ``src/waveDataFromSensorData_pluto.jl`` is a Pluto notebook that generates wave data of an irregular wave from measured sensor data. 
* ``examples/`` contains an example execution of the code (``exampleProgram.jl``)
* ``examples/results/`` contains the the scripts to compute the data (``computationScripts/``) used to plot all the figures (``plottingScripts/``) of the thesis which were not manually drawn.
* The [sensor data](https://doi.org/10.15480/882.9601) in form of a folder ``Data_Sensors/`` with subfolders ``With_Bathymetry`` and ``Without_Bathymetry`` needs to be copied into the ``examples/`` folder if new wave data is to be computed or the sensor data is to be plotted.

<!-- [![Build Status](https://github.com/ctech3000/InhomogeneousBathymetryEuler.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/ctech3000/InhomogeneousBathymetryEuler.jl/actions/workflows/CI.yml?query=branch%3Amaster) -->
