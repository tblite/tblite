# Light-weight tight-binding framework

[![License](https://img.shields.io/github/license/awvwgk/tblite)](https://github.com/awvwgk/tblite/blob/master/COPYING.LESSER)
[![Build Status](https://github.com/awvwgk/tblite/workflows/CI/badge.svg)](https://github.com/awvwgk/tblite/actions)
[![Documentation Status](https://readthedocs.org/projects/tblite/badge/?version=latest)](https://tblite.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/awvwgk/tblite/branch/main/graph/badge.svg?token=JXIE6myqNH)](https://codecov.io/gh/awvwgk/tblite)

This project is an effort to create a library implementation of the
extended tight binding (xTB) Hamiltonian which can be shared between
[``xtb``](https://github.com/grimme-lab/xtb) and
[``dftb+``](https://github.com/dftbplus/dftbplus).
The current state of this project should be considered as *highly experimental*.

Goals of this project are

- create a high-level interface to the extended tight binding methods
- allow low-level access to the components forming the actual energy expression
- provide a framework to handle and manipulate parametrisation data

Explicit non-goals are

- provide functionality beyond singlepoint calculations in this library
  (like geometry optimization or molecular dynamics)


## Installation

### Building from source

To compile this version of *tblite* the following programs are needed
(the number in parentheses specifies the tested versions).

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.55 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer
- a LAPACK / BLAS provider, like MKL or OpenBLAS

Optional dependencies are
- asciidoctor to build the manual page
- FORD to build the developer documentation
- C compiler to test the C-API and compile the Python extension module
- Python 3.6 or newer with the CFFI package installed to build the Python API

Setup a build with

```sh
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable.
To compile and run the projects testsuite use

```sh
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```sh
meson configure _build --prefix=/path/to/install
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.


## Usage

This project provides multiple entry points for different usage scenarios.
The simplest way to check out this project is by using the command line driver.


### Command line interface

The ``tblite`` runner executable provides full access to the implemented Hamiltonians.
You can run a single point calculation by providing a geometry input with

```
tblite --method gfn2 coord
```

For more details checkout the [``tblite(1)``](man/tblite.1.adoc) man page.


## Documentation

User and developer documentation is available [here](https://tblite.readthedocs.io).


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
Lesser GNU General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Lesser GNU General Public license, shall be licensed as above, without any
additional terms or conditions.
