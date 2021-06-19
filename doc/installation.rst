.. _install:

Installation
============

This project is currently in a highly experimental stage, therefore there are no stable distributions to offer yet.


Building from source
--------------------

This library depends on several Fortran modules to provide the desired functionality

- `mctc-lib`_: Modular computation tool chain library
- `dftd4`_: Reference implementation of the generally applicable charge-dependent London-dispersion correction, DFT-D4
- `s-dftd3`_: Reimplementation of the DFT-D3 dispersion correction
- `mstore`_: Molecular structure store (testing only)

.. _dftd4: https://github.com/dftd4/dftd4
.. _s-dftd3: https://github.com/awvwgk/simple-dftd3
.. _multicharge: https://github.com/grimme-lab/multicharge
.. _mctc-lib: https://github.com/grimme-lab/mctc-lib
.. _mstore: https://github.com/grimme-lab/mstore

.. _meson: https://mesonbuild.com
.. _ninja: https://ninja-build.org
.. _asciidoctor: https://asciidoctor.org
.. _cmake: https://cmake.org
.. _fpm: https://github.com/fortran-lang/fpm


Meson based build
~~~~~~~~~~~~~~~~~

The primary build system of this project is `meson`_.
For the full feature-complete build it is highly recommended to perform the build and development with meson.
To setup a build the following software is required

- A Fortran 2008 compliant compiler, like GCC Fortran and Intel Fortran classic
- `meson`_, version 0.55 or newer, better version 0.57.2 for the improved Fortran support
- `ninja`_, version 1.7 or newer, better version 1.10 for the dynamic dependency support
- a linear algebra backend, like MKL or OpenBLAS

Optional dependencies are

- `asciidoctor`_, to build the manual pages
- A C compiler to test the C API

To setup a new build run

.. code:: text

   meson setup _build --prefix=$HOME/.local

The Fortran and C compiler can be selected with the ``FC`` and ``CC`` environment variable, respectively.
For Intel Fortran oneAPI builds with MKL backend the ``-Dfortran_link_args=-qopenmp`` option has to be added.
To produce statically linked binaries set ``--default-library=static`` and add ``-Dfortran_link_args=-static`` as well.
The installation location is selected using the ``--prefix`` option.
The required Fortran modules will be fetched automatically from the upstream repositories and checked out in the *subprojects* directory.

To compile the project run

.. code:: text

   meson compile -C _build

Verify that the resulting projects is working correctly by running the testsuite with

.. code:: text

   meson test -C _build

In case meson is spawning to main concurrent test jobs limit the number of processes with the ``--num-processes`` option when running the test command.
By default the whole library and its subprojects are tested, to limit the testing to the project itself add ``--suite tblite`` as option.

Finally, you can make the project available by installation with

.. code:: text

   meson install -C _build


CMake based build
~~~~~~~~~~~~~~~~~

This project also provides support for `CMake`_ to give projects using it as build system an easier way to interface.
The CMake build files usually do not provide a feature-complete build, but contributions are more than welcome.
To setup a build the following software is required

- A Fortran 2008 compliant compiler, like GCC Fortran and Intel Fortran classic
- `cmake`_, version 3.14 or newer
- `ninja`_, version 1.10 or newer
- a linear algebra backend, like MKL or OpenBLAS

Configure a new build with

.. code:: text

   cmake -B _build -G Ninja -DCMAKE_INSTALL_PREFIX=$HOME/.local

You can set the Fortran compiler in the ``FC`` environment variable.
The installation location can be selected with the ``CMAKE_INSTALL_PREFIX``, GNU install directories are supported by default.
CMake will automatically fetch the required Fortran modules, you can provide specific version in the *subprojects* directory which will be used instead.

To run a build use

.. code:: text

   cmake --build _build

After a successful build make sure the testsuite passes

.. code:: text

   pushd _build && ctest --output-on-failure && popd

To make the project available install it with

.. code:: text

   cmake --install _build


Fpm based build
~~~~~~~~~~~~~~~

This projects supports building with the Fortran package manager (`fpm`_).
Create a new build by running

.. code:: text

   fpm build

You can adjust the Fortran compiler with the ``--compiler`` option and select the compilation profile with ``--profile release``.
To test the resulting build run the testsuite with

.. code:: text

   fpm test

The command line driver can be directly used from fpm wih

.. code:: text

   fpm run --profile release -- --help

To make the installation accessible install the project with

.. code:: text

   fpm install --profile release --prefix $HOME/.local
