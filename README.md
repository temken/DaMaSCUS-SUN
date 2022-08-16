[![Build & Test Status](https://github.com/temken/DaMaSCUS-SUN/workflows/Build%20&%20Tests/badge.svg)](https://github.com/temken/DaMaSCUS-SUN/actions)
[![codecov](https://codecov.io/gh/temken/DaMaSCUS-SUN/branch/main/graph/badge.svg)](https://codecov.io/gh/temken/DaMaSCUS-SUN)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# DaMaSCUS-SUN

<a href="https://ascl.net/2102.018"><img src="https://img.shields.io/badge/ascl-2102.018-blue.svg?colorB=262255" alt="ascl:2102.018" /></a>
[![DOI](https://zenodo.org/badge/263334878.svg)](https://zenodo.org/badge/latestdoi/263334878)
[![arXiv](https://img.shields.io/badge/arXiv-2102.12483-B31B1B.svg)](https://arxiv.org/abs/2102.12483)

Dark Matter Simulation Code for Underground Scatterings - Sun Edition

<img width="350" src="https://user-images.githubusercontent.com/29034913/108851182-73c1b500-75e4-11eb-94fd-11b93a3af0ae.png">

DaMaSCUS-SUN is a Monte Carlo tool simulating the process of solar reflection of dark matter (DM) particles as described in detail in [this publication](https://arxiv.org/abs/2102.12483).

## General Notes

- Solar reflection is the process of a DM particle from the galactic halo falling into the gravitational well of the Sun, where it can scatter on hot nuclei or electrons. If the DM particle scatters and escapes the Sun we consider it as "reflected".
- With DaMaSCUS-SUN we can describe this process by simulating individual DM particle's trajectories inside the Sun and in the solar system.
- Based on the Monte Carlo simulations, we can estimate the properties of the dark particle flux ejected from the Sun and derive e.g. exclusion limits.
- DaMaSCUS-SUN is written in C++ and built with CMake.
- The code is fully parallelized with MPI and can run on HPC clusters.

For more physics details, we refer to the [paper](https://arxiv.org/abs/2102.12483).

<details><summary>Repository content</summary>
<p>

The included folders are:

- *bin/*: This folder contains the executable after successful installation together with the configuration files.
- *data/*: Contains additional data necessary for the simulations, e.g. the solar model tables.
- *external/*: This folder will only be created and filled during the build with CMake and will contain the [obscura](https://github.com/temken/obscura) library necessary for all direct detection computations.
- *include/*: All header files of DaMaSCUS-SUN can be found here.
- *results/*: Each run of DaMaSCUS-SUN generates result files in a dedicated sub-folder named after the run's simulation ID string, which is specified in the configuration file.
- *src/*: Here you find the source code of DaMaSCUS-SUN.
- *tests/*: All code and executable files of the unit tests are stored here.

</p>
</details>

## Getting started

<details><summary>1. Dependencies</summary>
<p>

Before we can install DaMaSCUS-SUN, we need to make sure that a few dependencies are taken care of.

- [CMake](https://cmake.org/): DaMaSCUS-SUN as well as the libraries libphysica and obscura are built with CMake.
- [boost](https://www.boost.org/): For numerical integration (used by libphysica).
- [libconfig](https://github.com/hyperrealm/libconfig): For the configuration files, DaMaSCUS-SUN uses the libconfig library (required version at least 1.6). 
- [libphysica](https://github.com/temken/libphysica): Automatically downloaded to */external/obscura/external/*, compiled, and linked by CMake.
- [obscura](https://github.com/temken/obscura): Automatically downloaded to */external/*, compiled, and linked by CMake.
- [open MPI](https://www.open-mpi.org/): For the parallelization DaMaSCUS-SUN uses the open Message Passing Interface (MPI). Open MPI can be installed on Macs using [homebrew](https://brew.sh/):


<details><summary>Installation of boost</summary>
<p>

```
>brew install boost
```

or alternatively with APT:

```
>sudo apt-get install libboost-all-dev
```

</p>
</details>

<details><summary>Installation of libconfig</summary>
<p>
It is no-longer strictly necessary to install libconfig. During the build of *libphysica*, CMake will download and build *libconfig* in libphysica/external/ if it cannot find an installation. It still makes sense to install it in most cases.

On Macs, *libconfig* can be on installed using [homebrew](https://brew.sh/)

```
>brew install libconfig
```

or using APT on Linux machines

```
>sudo apt-get update -y
>sudo apt-get install -y libconfig++-dev
```

Alternatively, it can be built from the source files via

```
>wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
>tar -xvzf libconfig-1.7.2.tar.gz
>pushd libconfig-1.7.2
>./configure
>make
>sudo make install
>popd
```

</p>
</details>

<details><summary>Installation of open MPI</summary>
<p>

```
>brew install open-mpi
```

or alternatively with APT:

```
>sudo apt-get install openmpi-bin libopenmpi-dev
```

</p>
</details>

</p>
</details>

<details><summary>2. Downlad & Installation</summary>
<p>
The DaMaSCUS-SUN source code can be downloaded by cloning this git repository:

```
>git clone https://github.com/temken/DaMaSCUS-SUN.git 
>cd DaMaSCUS-SUN
```

The code is compiled and the executable is created using CMake.

```
>cmake -E make_directory build
>cd build
>cmake -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF ..
>cmake --build . --config Release
>cmake --install .
```

If everything worked well, there should be the executable *DaMaSCUS-SUN* in the */bin/* folder.

</p>
</details>

<details><summary>3. Usage</summary>
<p>
Once DaMaSCUS-SUN is installed, it can run by running the following command from the */bin/* folder:

```
>mpirun -n N DaMaSCUS-SUN config.cfg
```

Here *N* must be specified to the number of MPI processes. 

```
>mpirun -n N DaMaSCUS-SUN config.cfg                                                                                                                                                                                                                        1 ↵  
[Started on Wed Feb 24 13:19:23 2021]
DaMaSCUS-SUN-0.1.0	git:main/757a821
  ___       __  __      ___  ___ _   _ ___     ___ _   _ _  _
 |   \ __ _|  \/  |__ _/ __|/ __| | | / __|___/ __| | | | \| |
 | |) / _` | |\/| / _` \__ \ (__| |_| \__ \___\__ \ |_| | .` |
 |___/\__,_|_|  |_\__,_|___/\___|\___/|___/   |___/\___/|_|\_|
                               developed by Timon Emken (2020)

MPI processes:	N

##############################################################
Summary of obscura configuration

Config file:	config.cfg
ID:		identifier

----------------------------------------
DM particle summary:
	Mass:			1 MeV
	Spin:			0.5
	Low mass:		[x]
...
```

But before running DaMaSCUS-SUN, the user has to specify the scenario to be simulated by DaMaSCUS-SUN by setting a number of input parameters in the configuration file *config.cfg*.

</p>
</details>


<details><summary>4. The configuration file</summary>
<p>

In the configuration file, the user decides what scenario of solar reflection to simulate with DaMaSCUS-SUN.
A number of parameters need to be specified and they are described here.


| Note, that the handling of the input files by libconfig is sensitive to the data type. For example, a DM mass of 1 GeV has to be set as "1.0", **not** "1". |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------- |

1. First of all, it makes sense to give the simulation run a unique ID. All results will be saved under */results/identifier/*.

```
//DaMaSCUS-SUN - Configuration File

//ID
	ID		=	"identifier";
```

2. Next we need to decide if we run a single parameter point, or scan a grid of parameters to find detection exclusion limits. In either case, we need to specify the minimal number of data points to be generated by DaMaSCUS-SUN.

```
//Run mode
	run_mode = "Parameter point";	//Options: "Parameter scan" or "Parameter point"

	sample_size 		=	100;
	interpolation_points	=	1000;	//The scattering rate is interpolated on a NxN grid to speed up the simulations.
						//Recommended value: 1000
						//Set to 0 to run without interpolation.
```

3. If we run a parameter scan and compute exclusion limits, we need to specify the parameter grid.

```
//Options for "Parameter point"
	isoreflection_rings 		=	3;

// Options for "Parameter scan"
	compute_halo_constraints	= 	true;
	perform_full_scan		=	false;	//Full scan or STA contour tracing
	
	constraints_certainty		=	0.95;	//Certainty level
	
	constraints_mass_min		=	2.0e-6;	//in GeV										
	constraints_mass_max		=	1.0e-3;	//in GeV
	constraints_masses		=	5;

	cross_section_min 		=	1.0e-37;// in cm*cm
	cross_section_max 		=	1.0e-32;// in cm*cm	
	cross_sections			=	5;

```

4. The next block determines the DM particle properties. In the case of a "Parameter point" run, we need to set the DM mass and cross sections.

```
//Dark matter particle
	DM_mass		  		=	0.001;	// in GeV
	DM_spin		  		=	0.5;
	DM_fraction			=	1.0;	// the DM particle's fractional abundance (set to 1.0 for 100%)
	DM_light			=	true;	// Options: true or false. low mass mode

	DM_interaction			=	"SI";	// Options: "SI", "SD", or "DP"

	DM_isospin_conserved		=	true;   // only relevant for SI and SD
	DM_relative_couplings		=	(1.0, 0.0); //relation between proton (left) and neutron (right) couplings (only relevant if 'DM_isospin_conserved' is false.)
	DM_cross_section_nucleon	=	1.0e-80;    //in cm^2
	DM_cross_section_electron	=	1.0e-35;    //in cm^2 (only relevant for SI and SD)
	DM_form_factor			=	"Contact";	// Options: "Contact", "Electric-Dipole", "Long-Range", "General" (only relevant for SI and DP)
	DM_mediator_mass		=	0.0;	// in MeV (only relevant if 'DM_form_factor' is "General")

```

5. For the direct detection experiment, we can either choose one of the pre-defined experimental analyses or define an experiment ourselves.

```
//Dark matter detection experiment
	DD_experiment	=	"SENSEI@MINOS";	//Options for nuclear recoils: "Nuclear recoil", "DAMIC-2012", "XENON1T-2017", "CRESST-II","CRESST-III", "CRESST-surface"
 										//Options for electron recoils: "Semiconductor","protoSENSEI@MINOS","protoSENSEI@surface", "SENSEI@MINOS", "CDMS-HVeV", "Ionization", "XENON10-S2", "XENON100-S2", "XENON1T-S2", "DarkSide-50-S2"

	//Options for user-defined experiments ("Nuclear recoil", "Ionization", and "Semiconductor")
	//General
	DD_exposure 		=	300.0;	//in kg years
	DD_efficiency 		=	1.0;	//flat efficiency
	DD_observed_events 	=	0;	//observed signal events
	DD_expected_background 	=	0.0;	//expected background events

	//Specific options for "Nuclear recoil"
	DD_targets_nuclear	=	(
					(4.0, 8),
					(1.0, 20),
					(1.0, 74)
					);	// Nuclear targets defined by atom ratio/abundances and Z
	DD_threshold_nuclear	=	0.1;	//in keV
	DD_Emax_nuclear		=	40.0;	//in keV
	DD_energy_resolution	=	0.0;	//in keV
	
	//Specific options for Ionization and Semiconductor:
	DD_target_electron	=	"Xe";	//Options for Ionization: 	Xe, Ar
			    	//Options for Semiconductor:	Si, Ge
	DD_threshold_electron	=	4; //In number of electrons or electron hole pairs.
```

6. The initial conditons of our simulations are sampled from the DM halo model, which is determined by the following parameters.

```
//Dark matter distribution
	DM_distribution     =	"SHM";  //Options: "SHM"
	DM_local_density    =	0.4;	//in GeV / cm^3
	
	//Options for "SHM"
	SHM_v0		=	220.0;			//in km/sec
	SHM_vObserver	=	(11.1, 232.2, 7.3);	//in km/sec
	SHM_vEscape	=	544.0;			//in km/sec

```

</p>
</details>

## Version History

- 23.02.2021: Release of version 0.1.0

## Everything else

<details><summary>Citing DaMaSCUS-SUN</summary>
<p>

If you decide to use this code, please cite the latest archived version,

> Emken, T., 2021, Dark Matter Simulation Code for Underground Scatterings - Sun Edition (DaMaSCUS-SUN) [Code, v0.1.1], Astrophysics Source Code Library, record [[ascl:2102.018]](https://ascl.net/2102.018), [[DOI:10.5281/zenodo.5957388]](https://doi.org/10.5281/zenodo.5957388)

Bibtex entry:

```
@software{DaMaSCUSsun,
  author = {Emken, Timon},
  title = {{Dark Matter Simulation Code for Underground Scatterings - Sun Edition~(DaMaSCUS-SUN) [Code, v0.1.1]}},
  year         = {2021},
  publisher    = {Zenodo},
  version      = {v0.1.1},
  doi          = {DOI:10.5281/zenodo.5957388},
  url          = {https://doi.org/10.5281/zenodo.5957388},
  howpublished={Astrophysics Source Code Library record \href{https://ascl.net/2102.018}{[ascl:2102.018]}. The code can be found under \url{https://github.com/temken/damascus-sun}. Version 0.1.1 is archived as \href{https://doi.org/10.5281/zenodo.5957388}{DOI:10.5281/zenodo.5957388}}
}
```

As well as the original publications,

> Emken, T. , 2021,  **Solar reflection of light dark matter with heavy mediators**, [Phys.Rev.D 105 (2022) 6, 063020](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.105.063020), [[arXiv:2102.12483]](https://arxiv.org/abs/2102.12483).

Bibtex entry:

```
@article{Emken:2021lgc,
    author = "Emken, Timon",
    title = "{Solar reflection of light dark matter with heavy mediators}",
    eprint = "2102.12483",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1103/PhysRevD.105.063020",
    journal = "Phys. Rev. D",
    volume = "105",
    number = "6",
    pages = "063020",
    year = "2022"
}
```

</p>
</details>

<details><summary>Author & Contact</summary>
<p>

The author of DaMaSCUS-SUN is Timon Emken.

For questions, bug reports or other suggestions please contact Timon Emken ([timon.emken@fysik.su.se](mailto:timon.emken@fysik.su.se)).
</p>
</details>

<details><summary>License</summary>
<p>

This project is licensed under the MIT License - see the LICENSE file.

</p>
</details>