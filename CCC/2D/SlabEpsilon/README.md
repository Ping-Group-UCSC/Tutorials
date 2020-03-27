Calculating Epsilon(z)
===================================

Prerequisites:
------------------------------------
* [QE](https://www.quantum-espresso.org/)
* [JDFTx](http://jdftx.org/)
* [Python](https://www.python.org/)
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](https://matplotlib.org/)


Instructions:
-----------------------------------

*Calculating the dielectric tensor (QE)*

1. We will begin with computing the dielectric tensor within Quantum ESPRESSO using `ph.x` and `dynmat.x`. Enter the `QE` directory.

```bash
cd QE
```

2. This folder contains input for `pw.x`, `ph.x`, and `dynmat.x`: `pw.in`, `ph.in`, and `dynmat.in`, respectively. The job script can be used to run the calculation on kairay. However you choose, run each of the calculations in the appropriate order (1. `pw.x`, 2. `ph.x`, 3. `dynmat.x`).

3. Open the file `dynmat.out`. Under the line `Electronic dielectric permittivity tensor (F/m units)` is the dielectric tensor from electron contributions only (eps_inf). Follwoing this, under the line ` ... with zone-center polar mode contributions`, is the dielectric tensor with both electron and ion contribution (eps_0).


```bash
Electronic dielectric permittivity tensor (F/m units)
         1.828286   -0.000000    0.000000
        -0.000000    1.828285    0.000000
         0.000000    0.000000    1.159872

 ... with zone-center polar mode contributions
         2.223744   -0.000000    0.000000
        -0.000000    2.223744   -0.000000
         0.000000   -0.000000    1.172014
```


*Calculating the perpendicular dielectric function (JDFTx)*

4. We will now compute the dielectric function. Exit the directory `QE` and enter the directory `JDFTx`

```bash
cd ../JDFTx
```

4. This folder contains two input files `minus.in` and `plus.in`. Will run both calculations and the `JDFTx` code will use their response to an external electric field to calculate the dielectric function. However you choose, run both calculations.

5. Open the file `plus_dump.slabEpsilon`, this contains the dielectric function due to a field perpendicular to our slab (eps_perp(z)).


*Calculating the parallel dielectric function (python)*

6. We will now calculate the parallel component of the dielectric function. Implementing the methods detailed in [Phys. Rev. Materials 1, 071001(R)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.1.071001), the python script `slab_eps.py` can be used to calculate the parallel compoenent of the dielectric function from the perpendicular component and the dielectric matrix which we computed above. Run the python script to generate the file `epsilon_z.dat`. It will also plot the data and save to the file `epsilon_z.eps`.

```bash
./slab_eps.py
```


Author(s)
------------------------------------
Tyler J. Smart

