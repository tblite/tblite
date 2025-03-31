Running calculations from Python
================================

The *tblite* Python package allows to run extended tight binding (xTB) calculations directly in Python.
This tutorial demonstrates how to set up and run a single-point calculation using GFN2-xTB.


Installing the package
----------------------

To start create a new Python environment using the mamba package manager.
We specify the packages we want to install in our environment file:

.. code-block:: yaml
   :caption: environment.yml

   name: xtb
   channels:
   - conda-forge
   dependencies:
   - tblite-python
   - qcelemental
   - polars

Save the file as *environment.yml* and create the environment by running

.. code-block:: bash

   mamba env create -n xtb -f environment.yml
   mamba activate xtb

This will create a new environment called *xtb* and install all the necessary packages.
Make sure that *tblite* is available in your Python environment.
You can check this by opening a Python interpreter and importing the package

.. code-block:: python

   In [1]: import tblite.interface as tb
   
   In [2]: tb.library.get_version()
   Out[2]: (0, 4, 0)


First calculation
-----------------

We will run our first calculation directly from Python.
To specify the input we will directly initialize the calculation without reading any external files.
For this we specify our input as cartesian coordinates in Bohr and atomic numbers for the elements.

.. code-block:: python

   In [1]: import numpy as np

   In [2]: coordinates = np.array([
      ...:     [ 2.02799738646442,  0.09231312124713, -0.14310895950963],
      ...:     [ 4.75011007621000,  0.02373496014051, -0.14324124033844],
      ...:     [ 6.33434307654413,  2.07098865582721, -0.14235306905930],
      ...:     [ 8.72860718071825,  1.38002919517619, -0.14265542523943],
      ...:     [ 8.65318821103610, -1.19324866489847, -0.14231527453678],
      ...:     [ 6.23857175648671, -2.08353643730276, -0.14218299370797],
      ...:     [ 5.63266886875962, -4.69950321056008, -0.13940509630299],
      ...:     [ 3.44931709749015, -5.48092386085491, -0.14318454855466],
      ...:     [ 7.77508917214346, -6.24427872938674, -0.13107140408805],
      ...:     [10.30229550927022, -5.39739796609292, -0.13672168520430],
      ...:     [12.07410272485492, -6.91573621641911, -0.13666499342053],
      ...:     [10.70038521493902, -2.79078533715849, -0.14148379504141],
      ...:     [13.24597858727017, -1.76969072232377, -0.14218299370797],
      ...:     [ 7.40891694074004, -8.95905928176407, -0.11636933482904],
      ...:     [ 1.38702118184179,  2.05575746325296, -0.14178615122154],
      ...:     [ 1.34622199478497, -0.86356704498496,  1.55590600570783],
      ...:     [ 1.34624089204623, -0.86133716815647, -1.84340893849267],
      ...:     [ 5.65596919189118,  4.00172183859480, -0.14131371969009],
      ...:     [14.67430918222276, -3.26230980007732, -0.14344911021228],
      ...:     [13.50897177220290, -0.60815166181684,  1.54898960808727],
      ...:     [13.50780014200488, -0.60614855212345, -1.83214617078268],
      ...:     [ 5.41408424778406, -9.49239668625902, -0.11022772492007],
      ...:     [ 8.31919801555568, -9.74947502841788,  1.56539243085954],
      ...:     [ 8.31511620712388, -9.76854236502758, -1.79108242206824],
      ...: ])

   In [3]: elements = np.array([6,7,6,7,6,6,6,8,7,6,8,7,6,6,1,1,1,1,1,1,1,1,1,1])

   In [4]: import tblite.interface as tb

   In [5]: xtb = tb.Calculator("GFN2-xTB", elements, coordinates)

   In [6]: results = xtb.singlepoint()
   ------------------------------------------------------------
     cycle        total energy    energy error   density error
   ------------------------------------------------------------
         1     -41.75162462696  -4.2243950E+01   1.9479957E-01
         2     -42.11867876340  -3.6705414E-01   7.5972202E-02
         3     -42.14180557544  -2.3126812E-02   4.6343403E-02
         4     -42.14537345276  -3.5678773E-03   1.2550676E-02
         5     -42.14691416477  -1.5407120E-03   6.1305240E-03
         6     -42.14742063287  -5.0646811E-04   2.4358092E-03
         7     -42.14744770792  -2.7075047E-05   1.0726515E-03
         8     -42.14746243589  -1.4727977E-05   3.6117558E-04
         9     -42.14746301451  -5.7861259E-07   1.6698189E-04
        10     -42.14746312007  -1.0556375E-07   5.4777956E-05
        11     -42.14746315728  -3.7209709E-08   2.5307046E-05
        12     -42.14746315838  -1.0997425E-09   8.2333731E-06
   ------------------------------------------------------------

    total:                                   2.318 sec

   In [7]: results.get("energy")
   Out[7]: array(-42.14746316)


While it is possible to run xTB calculations this way directly in Python, it will become quickly cumbersome if we want to run many calculations at once.
Instead we want to read our geometry from an input file, for example an xyz geometry file.

.. code-block:: none
   :caption: caffeine.xyz

   24

   C            1.07317        0.04885       -0.07573
   N            2.51365        0.01256       -0.07580
   C            3.35199        1.09592       -0.07533
   N            4.61898        0.73028       -0.07549
   C            4.57907       -0.63144       -0.07531
   C            3.30131       -1.10256       -0.07524
   C            2.98068       -2.48687       -0.07377
   O            1.82530       -2.90038       -0.07577
   N            4.11440       -3.30433       -0.06936
   C            5.45174       -2.85618       -0.07235
   O            6.38934       -3.65965       -0.07232
   N            5.66240       -1.47682       -0.07487
   C            7.00947       -0.93648       -0.07524
   C            3.92063       -4.74093       -0.06158
   H            0.73398        1.08786       -0.07503
   H            0.71239       -0.45698        0.82335
   H            0.71240       -0.45580       -0.97549
   H            2.99301        2.11762       -0.07478
   H            7.76531       -1.72634       -0.07591
   H            7.14864       -0.32182        0.81969
   H            7.14802       -0.32076       -0.96953
   H            2.86501       -5.02316       -0.05833
   H            4.40233       -5.15920        0.82837
   H            4.40017       -5.16929       -0.94780

Instead of implementing our own xyz file reader we will be using the qcelemental package which already provides this functionality for us.
Fortunately, the qcelemental library does store the geometry already in Bohr and we do not need to convert the coordinates to input them in our xTB calculation.

.. code-block:: python

   In [1]: import qcelemental as qcel

   In [2]: molecule = qcel.models.Molecule.from_file("caffeine.xyz")

   In [3]: import tblite.interface as tb

   In [4]: xtb = tb.Calculator("GFN2-xTB", molecule.atomic_numbers, molecule.geometry)

   In [5]: results = xtb.singlepoint()
   ------------------------------------------------------------
     cycle        total energy    energy error   density error
   ------------------------------------------------------------
         1     -41.75162462526  -4.2243950E+01   1.9479957E-01
         2     -42.11867876267  -3.6705414E-01   7.5972202E-02
         3     -42.14180557495  -2.3126812E-02   4.6343403E-02
         4     -42.14537345231  -3.5678774E-03   1.2550676E-02
         5     -42.14691416431  -1.5407120E-03   6.1305240E-03
         6     -42.14742063242  -5.0646811E-04   2.4358092E-03
         7     -42.14744770747  -2.7075046E-05   1.0726515E-03
         8     -42.14746243545  -1.4727977E-05   3.6117558E-04
         9     -42.14746301406  -5.7861259E-07   1.6698189E-04
        10     -42.14746311962  -1.0556377E-07   5.4777956E-05
        11     -42.14746315683  -3.7209702E-08   2.5307046E-05
        12     -42.14746315793  -1.0997567E-09   8.2333732E-06
   ------------------------------------------------------------
   
    total:                                   2.067 sec
   
   In [6]: results.get("energy")
   Out[6]: array(-42.14746316)

We find the calculation for the caffeine molecule is run and the energy we obtained before is also found again.


Evaluating properties
---------------------

Now that we can evaluate energies we want to extend the evaluation to other properties with xTB.
Let's compute the vertical ionization potential for caffeine with xTB.
For computing this property we have two options, first we can get the ionization potential directly from our xTB wavefunction by using the energy of the highest occupied orbital.
We continue from our previous session and obtain the orbital energies and occupation numbers to find the highest occupied orbital.
Its energy approximates the negative vertical ionization potential.

.. code-block:: python

   In [7]: orbital_energies = results.get("orbital-energies")

   In [8]: orbital_occupations = results.get("orbital-occupations")

   In [9]: import numpy as np

   In [10]: homo_index = np.argmin(orbital_occupations) - 1

   In [11]: -orbital_energies[home_index] * qcel.constants.conversion_factor("hartree", "eV")
   Out[11]: np.float64(10.592677756988177)


Instead of approximating the ionization potential we can also compute it.
The energy of removing an electron can be expressed by the reaction

.. math::

   \text{Molecule} \leftarrow \text{Molecule}^{+} + e^{-}

This reaction energy is the negative ionization potential.
To compute this energy with xTB we update our calculator by setting the total charge to +1:

.. code-block:: python

   In [12]: xtb.update(charge=1)

   In [13]: results_ion = xtb.singlepoint()
   ------------------------------------------------------------
     cycle        total energy    energy error   density error
   ------------------------------------------------------------
         1     -41.26491540599  -4.1757240E+01   1.9757517E-01
         2     -41.55259512561  -2.8767972E-01   8.4795450E-02
         3     -41.62346372317  -7.0868598E-02   6.4944961E-02
         4     -41.66084952611  -3.7385803E-02   2.1386398E-02
         5     -41.65752591189   3.3236142E-03   1.3708257E-02
         6     -41.66153978913  -4.0138772E-03   8.1833393E-03
         7     -41.66258949427  -1.0497051E-03   4.0653874E-03
         8     -41.66278148206  -1.9198779E-04   2.3414702E-03
         9     -41.66288309037  -1.0160831E-04   9.2578279E-04
        10     -41.66289562556  -1.2535187E-05   4.9559562E-04
        11     -41.66289374167   1.8838873E-06   4.9598725E-04
        12     -41.66290207235  -8.3306780E-06   1.4008496E-04
        13     -41.66290260550  -5.3315508E-07   4.0745834E-05
        14     -41.66290259015   1.5358552E-08   4.3706476E-05
        15     -41.66290264433  -5.4182017E-08   5.4270862E-06
   ------------------------------------------------------------
   
    total:                                   3.178 sec

   In [14]: (results_ion.get("energy") - results.get("energy")) * qcel.constants.conversion_factor("hartree", "eV")
   Out[14]: np.float64(13.185563185488823) 


We do find quite a difference in the calculated value and the approximated one.
Before we can use the ionization potential computed by xTB we should however correct for the self-interaction error using an empirical determined shift of 4.846V.
This shift should be applied for all ionization potentials computed with xTB.

.. tip::

   Since xTB is a semiempirical method it makes some approximations which result in a strong self-interaction for a free electron.
   This value can be computed exactly from the xTB parameters or determined empirically.
   For a full derivation checkout Ref. :footcite:`neugebauer2020`.


Fukui indices from partial charges
----------------------------------

While the molecular ionization potential is a great descriptor for the whole molecule, xTB also provides many properties which are atom resolved
Computing the Fukui index provides a simple descriptor for chemical reactivity which we can compute from the partial charges according to the following equations: 

.. math::

   f_\text{A}^{(+)} = q_\text{A}^{(0)} - q_\text{A}^{(-1)} \quad
   f_\text{A}^{(-)} = q_\text{A}^{(+1)} - q_\text{A}^{(0)} \quad
   f_\text{A}^{(0)} = \frac12 \left(q_\text{A}^{(+1)} - q_\text{A}^{(-1)}\right)

where we have the three Fukui indices computed from the partial charges of the neutral (0), cationic (+) and anionic (-) system.
To perform this calculation with xTB we go back to our computation environment and update our molecule to a negative total charge:

.. code-block:: python

   In [15]: xtb.update(charge=-1)
   
   In [16]: results_neg = xtb.singlepoint()
   ------------------------------------------------------------
     cycle        total energy    energy error   density error
   ------------------------------------------------------------
         1     -41.90273954376  -4.2395064E+01   1.8739068E-01
         2     -42.14847705596  -2.4573751E-01   9.6117213E-02
         3     -42.22375508217  -7.5278026E-02   7.0987678E-02
         4     -42.30265756019  -7.8902478E-02   3.4694875E-02
         5     -42.30270404731  -4.6487120E-05   2.2615887E-02
         6     -42.31014040726  -7.4363600E-03   1.1765251E-02
         7     -42.31193645170  -1.7960444E-03   3.5463976E-03
         8     -42.31192828337   8.1683241E-06   1.3330928E-03
         9     -42.31195331590  -2.5032523E-05   5.6195325E-04
        10     -42.31195526941  -1.9535110E-06   2.4623254E-04
        11     -42.31195551932  -2.4991225E-07   1.1604168E-04
        12     -42.31195568210  -1.6278454E-07   3.6907929E-05
        13     -42.31195568751  -5.4046154E-09   1.4924180E-05
   ------------------------------------------------------------
   
    total:                                   2.845 sec

   In [17]: import polars as pl

   In [18]: fukui = pl.DataFrame({
       ...:     "element": molecule.symbols,
       ...:     "f(+)": results.get("charges") - results_neg.get("charges"),
       ...:     "f(-)": results_ion.get("charges") - results.get("charges"),
       ...:     "f(0)": (results_ion.get("charges") - results_neg.get("charges"))/2,
       ...: })

   In [19]: fukui
   Out[19]: 
   shape: (24, 4)
   ┌─────────┬───────────┬───────────┬──────────┐
   │ element ┆ f(+)      ┆ f(-)      ┆ f(0)     │
   │ ---     ┆ ---       ┆ ---       ┆ ---      │
   │ str     ┆ f64       ┆ f64       ┆ f64      │
   ╞═════════╪═══════════╪═══════════╪══════════╡
   │ C       ┆ -0.029669 ┆ -0.025591 ┆ -0.02763 │
   │ N       ┆ 0.070509  ┆ 0.060968  ┆ 0.065738 │
   │ C       ┆ 0.054677  ┆ 0.025012  ┆ 0.039845 │
   │ N       ┆ 0.0796    ┆ 0.068182  ┆ 0.073891 │
   │ C       ┆ 0.036333  ┆ 0.020819  ┆ 0.028576 │
   │ …       ┆ …         ┆ …         ┆ …        │
   │ H       ┆ 0.03337   ┆ 0.067995  ┆ 0.050683 │
   │ H       ┆ 0.033354  ┆ 0.067899  ┆ 0.050627 │
   │ H       ┆ 0.02604   ┆ 0.031893  ┆ 0.028967 │
   │ H       ┆ 0.038735  ┆ 0.040089  ┆ 0.039412 │
   │ H       ┆ 0.038748  ┆ 0.040102  ┆ 0.039425 │
   └─────────┴───────────┴───────────┴──────────┘



Literature
----------

.. footbibliography::