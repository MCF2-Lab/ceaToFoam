# Installation

## Requirements

| Component | Version | Used for |
| --- | --- | --- |
| Python | 3.9 or newer | everything |
| [Cantera](https://cantera.org) | 3.0 or newer | thermodynamic and transport data of the frozen mixture |
| [CEA_Wrap](https://github.com/civilwargeeky/CEA_Wrap) | 2.1 or newer | calling NASA CEA |
| NumPy, SciPy, PyYAML | any recent | numerics and input files |
| [OpenFOAM](https://www.openfoam.com) (openfoam.com) | tested with v2606 | running the generated case |
| Matplotlib | optional | the figure of `ceatofoam sweep --plot` |

The Python dependencies are installed by `pip`. OpenFOAM is needed only to run
the case, not to write it.

## Installing ceaToFoam

```sh
git clone https://github.com/MCF2-Lab/ceaToFoam.git
cd ceaToFoam
pip install .              # or: pip install '.[plot]' to include Matplotlib
```

A virtual environment keeps the dependencies separate from the system Python:

```sh
python3 -m venv ~/venvs/ceatofoam
. ~/venvs/ceatofoam/bin/activate
pip install '.[plot]'
```

## The CEA executable on Linux and macOS

CEA_Wrap ships the CEA executable for Windows only. On Linux and macOS,
compile the Fortran source that comes with it once (gfortran is enough) and
tell CEA_Wrap where the executable is:

```sh
ASSETS=$(python3 -c "import CEA_Wrap, os; print(os.path.join(os.path.dirname(CEA_Wrap.__file__), 'assets'))")
mkdir -p ~/.local/bin && cd "$ASSETS"
gfortran -O2 -std=legacy -o ~/.local/bin/Dan_CEA2 dan_cea2.f
export CEA_EXE_LOCATION=~/.local/bin/Dan_CEA2    # add to ~/.bashrc
```

Run the command in the environment where ceaToFoam is installed, so that
`python3` finds CEA_Wrap.

## Checking the installation

```sh
ceatofoam --version
ceatofoam design examples/rp1-lox-500psi.yaml
python3 -m unittest discover -s tests
```

The design command prints the chamber state and the design point of the
example engine (see [Usage](usage.md)). Tests that need the CEA executable or
OpenFOAM are skipped when these are not available; with both present, all
tests run, including `blockMesh` and `checkMesh` on a generated case.
