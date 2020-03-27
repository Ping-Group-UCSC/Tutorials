Charged Cell Correction of 2D systems
===================================

Description
------------------------------------
This repository provides a tutorial for calculating the charged cell correction of 2D systems. 
This involves first computing the dielectric profile `eps(z)` for the two-dimensional system.
Note! This tutorial assumes a user already has confidence in running QE calculations, as such all QE calculations (besides the calculation of the dielectric tensor) are already completed. Feel free to rerun the QE calculations with the given files if desired, but the focus here is on the steps needed to calculate CCC after one has calculated a few charged defect calculations in QE.

---

Prerequisites:
------------------------------------
* [QE](https://www.quantum-espresso.org/)
* [JDFTx](http://jdftx.org/)
* [Python](https://www.python.org/)
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](https://matplotlib.org/)

---

Outline:
------------------------------------
**Steps 1-6** Calculate the dielectric function eps(z) due to perpendicular and parallel fields: `epsilon_z.dat`.

**Steps 7-9** Calculate the total electrostatic potential of the pristine system: `prist.d_tot`.

**Steps 10-13** Calculate the charged cell correction (CCC) for the Ti defect in h-BN. 

---

Instructions:
------------------------------------

***Dielectric Profile***

0. Enter the directory `SlabEpsilon`.

```bash
cd SlabEpsilon
```

*Calculating the dielectric tensor (QE)*

1. We will begin by computing the dielectric tensor within Quantum ESPRESSO using `ph.x` and `dynmat.x`. Enter the `QE` directory.

```bash
cd QE
```

2. This folder contains input for `pw.x`, `ph.x`, and `dynmat.x`: `pw.in`, `ph.in`, and `dynmat.in`, respectively. The script `job` can be used to run each calculation on Kairay. However you choose, run each of the calculations in the appropriate order (1. `pw.x`, 2. `ph.x`, 3. `dynmat.x`).

```bash
sbatch job
```

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

4. This folder contains two input files `minus.in` and `plus.in`. We will run both calculations with the `JDFTx` code (1. `minus.in`, 2. `plus.in`). However you choose, run both calculations.

```bash
sbatch job
```

5. Open the file `plus_dump.slabEpsilon`, this contains the dielectric function due to a field perpendicular to our slab (eps_perp(z)). After you are done, exit the `JDFTx` directory.

```bash
cd ..
```

*Calculating the parallel dielectric function (python)*

6. We will now calculate the parallel component of the dielectric function. The python script `slab_eps.py` can be used to calculate the parallel component of the dielectric function from the perpendicular component and the dielectric matrix which we computed above. See methods detailed in [Phys. Rev. Materials 1, 071001(R)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.1.071001). Run the python script to generate the file `epsilon_z.dat`. It will also plot the data and save to the file `epsilon_z.eps`.

```bash
./slab_eps.py
```

***Charge Cell Correction (CCC)***

*Pristine JDFTx calculation*

7. Exit the directory `SlabEpsilon` and enter the `JDFTx` directory within the pristine folder `Prist`.

```bash
cd ../Prist/JDFTx
```

8. Before we calculate the charge cell correction it is necessary to complete a calculation of the pristine system. Within this directory we have the basic set up of a JDFTx calculation for pristine h-BN. Note that we need the cell dimensions of this calculation to match for our later calculations so we cannot simply do a unit cell calculation.  However you choose, run the jdftx calculation.

```bash
sbatch job
```

9. This calculation generates the important file `prist.d_tot`, corresponding to the total electrostatic potential of the pristine system. This file alongside `epsilon_z.dat` from step 6 are all we need from here on out to calculate the charge cell correction for defects in h-BN. Note that if we change the cell dimensions (supercell size) or if we change the cutoff energy, we will have to repeat steps 1-8 and create new `prist.d_tot` and `epsilon_z.dat`.

*Charge cell correction for the Ti defect in h-BN*

10. Exit the directory `Prist/JDFTx` and enter the directory `Ti`.

```bash
cd ../../Ti
```

11. In this directory is four different charged calculations of the Ti defect (`Q+1`, `Q+2`, `Q-1`, `Q-2`) alongside the neutral defect calculation `Q0`. Directly within each folder is a completed QE calculation `relax.in` and `relax.out`. Using the command `grep tot_charge */relax.in` we can see that each calculation was given a different charge.

```bash
grep tot_charge */relax.in
```

12. For the cases where `tot_charge != 0`, there is a subfolder `CCC` which has the necessary input to run a JDFTx calculation to calculate the charge cell correction (CCC) using the total electrostatic potential of the pristine system `prist.d_tot` (computed in step 8) and the dielectric profile `epsilon_z.dat` (computed in step 6). If you open the file `Q+1/CCC/main.in` for example you will find the following three lines:

```bash
elec-initial-charge -1 # -1 electrons -> q= +1
charged-defect 0.490470496 0.509445064 0.579405106 -1 1 # location of Ti defect with charge -1
charged-defect-correction Slab 001 /export/data/share/tjsmart/h-BN/2019-Various/TiBN/QE/S07/Prist/JDFTx/prist.d_tot epsilon_z.dat  8. 1.
```

13. However you choose, run each of the calculations under the `Q*/CCC` folders.

```bash
hdir=$PWD
for d in $(find . -name 'CCC'); do
        cd $d
        sbatch job
        cd $hdir
done
```

***Post Processing and Plotting***

14. *Coming soon!*


Author(s)
------------------------------------
Tyler J. Smart

