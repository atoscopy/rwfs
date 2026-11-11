# README

# RWFS
## An Open-Source Python Package for Bound-State Solutions of the 
## Radial Dirac Equation in Arbitrary Central-Field Potentials

**Authors**
 M. Z. Mansour, Yosr E. E.-D. Gamal, S. Fritzsche,Th. M. El Sherbini
 
E-mail: **mmansour@niles.edu.eg**

----------------------------------------------------------------------

## Overview

RWFS is an open-source Python package for solving the radial Dirac
equation for bound electrons in arbitrary central-field potentials.
The program computes relativistic bound-state orbitals, orbital
energies, radial expectation values, and normalized radial
wavefunctions for atoms and ions.

Several analytical model potentials are available, including

- Thomas–Fermi,
- Coulomb,
- Yukawa,
- Multi-Yukawa,
- Screened,
- Screened-hydrogenic,
- finite-nuclear-size models: uniformly charged sphere and Fermi
  charge distribution,

as well as user-supplied numerical potentials read directly from an
external file (nine central-field potentials in total).

RWFS is restricted to bound states of electrons in attractive central
potentials (energies -c^2 < E < 0, measured from the rest-mass threshold).
The negative-energy region of the Dirac spectrum (E <= -2c^2) and the
positive-energy continuum are not treated.

The calculated radial wavefunctions can be written in ASCII format and
converted into the binary format required by the GRASP package for
subsequent Multiconfiguration Dirac–Hartree–Fock (MCDHF),
Relativistic Configuration Interaction (RCI), and related atomic
structure calculations.

----------------------------------------------------------------------

## Directory Structure

The RWFS distribution contains the following directories

RWFS/
│
├──> src/
├──> examples/
├──> tools/
├──> LICENSE
└──> README.md

----------------------------------------------------------------------

## src/

The 'src' directory contains the RWFS source code.

Main files include

- 'rwfs.py'
  Main program.

- 'energy_loop.py'
  Radial Dirac solver and energy iteration.

- 'grid.py'
  Generation of the numerical radial grid.

- 'potential.py'
  Construction of the selected central-field potential, including the
  finite-nuclear-size models.

- 'spline.py'
  Cubic-spline interpolation and numerical integration.

- 'mathematic.py'
  Numerical linear-algebra routines.

- 'constants.py'
  Physical constants and numerical parameters (including the Bohr radius
  in fm and the default Fermi skin thickness).

Running

python rwfs.py


starts either an interactive calculation or a batch calculation,
depending on the supplied input.

----------------------------------------------------------------------

## examples/

The 'examples' directory contains complete working examples illustrating
the use of RWFS.

The examples demonstrate both

- interactive execution, and
- fully automated batch calculations.

Typical example directories include

example_1/
example_2/
example_3/
example_4/
example_6/   (finite-nuclear-size effects, H to Db104+)

Each example contains the files required to reproduce the calculation,
such as

- input files,
- output files,
- summary files,
- radial wavefunctions,
- GRASP-compatible wavefunctions,
- generated potentials.

These examples are intended as templates for new calculations and
provide reproducible demonstrations of the package capabilities.

----------------------------------------------------------------------

## Central-field potentials

The potential is selected by number at the prompt (or in the batch input):

| Option | Potential |
|---|----------------|
| 1 | Thomas-Fermi |
| 2 | Coulomb |
| 3 | Yukawa |
| 4 | Multi-Yukawa |
| 5 | Screened Coulomb |
| 6 | Screened hydrogenic |
| 7 | User-defined potential read from a two-column (r, V(r)) file |
| 8 | Uniformly charged sphere (finite nucleus) |
| 9 | Fermi charge distribution (finite nucleus) |

### Finite nuclear size (options 8 and 9)

Options 8 and 9 give the potential of an extended nuclear charge for a given
nuclear rms charge radius R_rms. Only the bare nuclear potential is
produced (no electron screening), as for the Coulomb option, so they are
meant for hydrogen-like ions or as the nuclear part of a model potential.

- **Option 8, uniformly charged sphere:** R = sqrt(5/3) R_rms.
- **Option 9, Fermi distribution:** rho(r) = rho_0 / (1 + exp[(r - c)/a]),
  with a = t/(4 ln 3) and skin thickness t = 2.3 fm (default); the
  half-density radius c is fixed by R_rms. This requires
  R_rms > 1.95 fm for t = 2.3 fm.

----------------------------------------------------------------------

## tools/

The 'tools' directory contains auxiliary utilities for post-processing
RWFS calculations.

Currently the directory includes

readrwf_grasp.f

This Fortran utility converts the ASCII radial wavefunctions generated
by RWFS into the binary wavefunction format required by the GRASP
package.

The resulting binary file can be used directly as input for GRASP
MCDHF, RCI, and other relativistic atomic-structure calculations.

Detailed compilation and usage instructions are provided in the
accompanying README for 'readrwf_grasp.f'.

----------------------------------------------------------------------

## Typical Workflow

The normal workflow is

Run RWFS
   |
   v
Calculate bound orbitals
   |
   v
ASCII radial wavefunctions (*.rwf)
   |
   v
readrwf_grasp.f
   |
   v
Binary GRASP wavefunctions (*.w)
   |
   v
GRASP MCDHF / RCI calculations

----------------------------------------------------------------------

## Requirements

- Python 3.0 or later

- GNU Fortran (optional, only for compiling 'readrwf_grasp.f')

## Installation and Usage

The package consists of multiple interdependent Python modules. To run the program:

1. Ensure that Python 3 is installed.
2. Open a terminal.
3. Execute the main script by specifying its path:

python /path/to/src/rwfs.py

Replace '/path/to/src/' with the location of the 'src' directory on your system.
