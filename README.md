<!--- BOF -->
## CORPAD

This is a C++ code intended for numerical simulations that help investigate the spatial diffusion of relativistic charged particles.
These "test" particles are propagated in a medium with non-zero magnetic field **Bo**, which has a modeled turbulent component with zero mean.
The collective behaviour of these test particles allow us to determine values of the mean free paths (i.e. diffusion coefficients), in parallel and perpendicular direction to **Bo**, as a function of the turbulence properties of the magnetic medium.


---
### Compilation steps

This should compile and build an executable called `CRs.diff.x`.
```bash
cd src/
make default
# this uilds a parallel code, so it nedds MPI libraries
# NOTE: this was tested with an C++11 compiler.
```

In order to run a simulation, you must prepare input files `inp_turb`, `inp_ori` and `inp_gral`, where each file contains:
* `inp_turb`: parameters for the modeled magnetic turbulence
* `inp_ori` : the initial orientations of the particles
* `inp_gral`: other parameters related to: the accuracy of the simulation, the number of turbulence realizations of the model, and monitoring of the scattering of the particles.
For a test run suite, check [this](scripts/test.sh) bash script, and the input files referenced therein.


---
_Documentation in progress..._


<!--- EOF -->
