<!--- BOF -->
## CORPAD

This is a C++ code intended for numerical simulations that help investigate the spatial diffusion of relativistic charged particles (i.e. Cosmic Rays in out context).

These "test" particles are propagated in a medium with non-zero magnetic field **Bo**, which has a modeled turbulent component with zero mean.
The collective behaviour of these test particles allow us to determine values of the mean free paths (i.e. diffusion coefficients), in parallel and perpendicular direction to **Bo**, as a function of the turbulence properties of the magnetic medium.


---
### Particle code

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


---
### Evaluation of the turbulent model
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
## Pre-processing

This stage is optional. 
It's only necessary if you want to automatically generate the `.in` input files using empirical properties of the solar wind that are dependent on the helioradius `ro`, for any given energy of the particles.

The script [gen_input.py](scripts/gen_input.py) must be manually modified and then run it, so the `.in` files are generated in the [inputs](inputs) directory.

The parameters available for modification in the script includes: energy of the particles, background magnetic field, heliodistance, turbulence properties, total number of particles, maximum simulation time, resolution on the histograms of particle scatterings, and many more.



---
## Post-processing

The simulation runs produce raw data of the individual behaviour of the particles.
In order to obtain the diffusion coefficients of the collective behaviour of the particles, you can use the Python scripts in [this](scripts/analyze_err/atol) directory; also the check the other [README](scripts/analyze_err/atol/README.md).
Usually, you'll mostly use [massive.k.vs.t.py](scripts/analyze_err/atol/massive.k.vs.t.py), which determines the mean free paths (parallel and perpendicular to **Bo**) as a function of time.

Check all the options with `./massive.k.vs.t.py -- -h`.
Note that the script already has default values for the seeds and `t_decr`, and it's often not necessary to specify different values for these in the command line:
```bash
$ ./calc_k_vs_t.py -- -h
usage: massive.k.vs.t.py [-h] [-di DIR_SRC] [-do DIR_DST] [-df DIR_FIG]
                         [-sm_pe SEED_M_PERP] [-sb_pe SEED_B_PERP]
                         [-sm_pa SEED_M_PARA] [-sb_pa SEED_B_PARA]
                         [-td T_DECR]

optional arguments:
  -h, --help            show this help message and exit
  -di DIR_SRC, --dir_src DIR_SRC
                        source dir, relative to repo root
  -do DIR_DST, --dir_dst DIR_DST
                        output dir, relative to repo root
  -df DIR_FIG, --dir_fig DIR_FIG
                        output dir (relative to repo root)
  -sm_pe SEED_M_PERP, --seed_m_perp SEED_M_PERP
                        seed value 'm' for hyperbola fit-function
  -sb_pe SEED_B_PERP, --seed_b_perp SEED_B_PERP
                        seed value 'b' for hyperbola fit-function
  -sm_pa SEED_M_PARA, --seed_m_para SEED_M_PARA
                        seed value 'm' for hyperbola fit-function
  -sb_pa SEED_B_PARA, --seed_b_para SEED_B_PARA
                        seed value 'b' for hyperbola fit-function
  -td T_DECR, --t_decr T_DECR
                        value 't_decr' for hyperbola fit-function. It gives
                        the low threshold value in time, from which the data
                        is fitted.
```



---
_Documentation in progress..._


<!--- EOF -->
