<!--- BOF -->
## CORPAD

This is a C++ code intended for numerical simulations that help investigate the spatial diffusion of relativistic charged particles.
These "test" particles are propagated in a medium with non-zero magnetic field **Bo**, which has a modeled turbulent component with zero mean.
The collective behaviour of these test particles allow us to determine values of the mean free paths (i.e. diffusion coefficients), in parallel and perpendicular direction to **Bo**, as a function of the turbulence properties of the magnetic medium.


---
### Code compilations

For the code that performs the collective propagation of charged particles:
```bash
cd src/
make default
# this builds a parallel code, so it nedds MPI libraries
# NOTE: this was tested with an C++11 compiler.
```
This should compile and build an executable called `CRs.diff.x` in the root directory of this repository.

In order to run a simulation, you must prepare input files `inp_turb`, `inp_ori` and `inp_gral`, where each file contains:
* `inp_turb`: parameters for the modeled magnetic turbulence
* `inp_ori` : the initial orientations of the particles
* `inp_gral`: other parameters related to: the accuracy of the simulation, the number of turbulence realizations of the model, and monitoring of the scattering of the particles.
For a test run suite, check [this](scripts/test.sh) bash script, and the input files referenced therein.


We can also evaluate the magnetic turbulence model (isolated), which uses the same source files as the particle code above.
To compile, just:
```bash
cd src_Bmodel/
make cython
```
That builds a python interface for the C++ code that implements the turbulence model.
Then the model can be evaluated by importing it in any Python script as:
```python
import src_Bmodel.Bmodel
```
For a functional test, see this [test.py](src_Bmodel/test.py) script.




---
_Documentation in progress..._


<!--- EOF -->
