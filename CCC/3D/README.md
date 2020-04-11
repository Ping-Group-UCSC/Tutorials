Charged Cell Correction of 3D systems
===================================

Description
------------------------------------
This repository provides a tutorial for calculating the charged cell correction of 3D systems.

---

Prerequisites:
------------------------------------
* [QE](https://www.quantum-espresso.org/)
* [JDFTx](http://jdftx.org/)
---

Instructions:
------------------------------------

---

*Single point SCF in JDFTx (LaFeO3)*

1. Enter the directory `LaFeO3/Prist` to start things off by running a pristine calculation of LaFeO3 in JDFTx.

```bash
cd LaFeO3/Prist
sbatch job
```

2. From this calculation we generate the total electrostatic potential file for the pristine system `prist.d_tot`. This file will be used in subsequent calculations. Note if we change anything about the FFT grid of the system (i.e. change cutoff, lattice parameter, or directly set a different FFT) then we will need to redo step 1 and recreate a new `prist.d_tot`.

3. Enter the directory `../K` where we have calculations of the K doped LaFeO3 system in three different charge states `Q0`, `Q-1`, and `Q-2`. Enter each directory and run the calculation. Note that the calculations `Q-1` and `Q-2` depend on step 1 being completed, so if you submit the job scripts make sure to set the `--dependency` flag appropriately within the job script file.

```bash
cd ../K
cd Q0 ; sbatch job
cd ../Q-1 ; sbatch job
cd ../Q-2 ; sbatch job
```

4. Open the input files in `Q-1` and `Q-2` to the electronic charge is set with the command `elec-initial-charge`, we determine the placement and again confirm the charge of the defect with the command `charged-defect` and finally we enable charge defect correciton with the command `charged-defect-correction` wherein we specify the `prist.d_tot` file we generated in step 1.

```bash
vim -p Q-*/main.in
```

look for these lines in in Q-1/main.in:

```
    elec-initial-charge +1 # +1 electrons -> q= -1
    charged-defect 0.498148191 0.237282612 0.374922180 +1 1
    charged-defect-correction /export/data/share/tjsmart/LaFeO3/Prist/JDFTx/prist.d_tot  29 8. 1.
```

look for these lines in in Q-2/main.in:

```
    elec-initial-charge +2 # +2 electrons -> q= -2
    charged-defect 0.499103562 0.236896782 0.374996811 +2 1
    charged-defect-correction /export/data/share/tjsmart/LaFeO3/Prist/JDFTx/prist.d_tot  29 8. 1.
```

5. Check the charge cell correction you obtained and compare with the below.

```bash
grep Net Q-*/main.out
```
Produces the output:

```
    Q-1/main.out:	Net correction: 0.00622852 (= EmodelIsolated - EmodelPeriodic + q deltaV)
    Q-2/main.out:	Net correction: 0.01488844 (= EmodelIsolated - EmodelPeriodic + q deltaV)
```

<br>

---

*Fixed charge density from QE in JDFTx (Fe2O3)*

1. The procedure here is essentially identical to above with a few added steps. The difference here is that we will not be rerunning any unecessary calculations in JDFTx. We will only simply be looking to get the charged defect correction value and we will import the charge density from QE. In order to convert the charge density we will use the script [conv.py](https://gitlab.com/ping-group-ucsc/qe-jdftx-convert). On kairay you can find this script under the following directory:

```bash
/export/data/share/tjsmart/Programs/Ping-Group/qe-jdftx-convert/conv.py
```

2. This calculation requires doing things in a slightly different order in order to ensure that the FFT grids match between all of our calculations (note: that while we set the same ecut in QE and JDFTx we do not always get the same fft grid!). First enter the directory `Fe2O3/Sn-Cluster`. We will run the script `run_conv.sh` in order to convert the charge density of the charge calculations in `Q+1` and `Q+2`. Enter both directories and run this script.

```bash
cd Q+1
../run_conv.sh
cd ../Q+2
../run_conv.sh
```

You should see the following output printed when the script is ran:

```
    Charge density dimension [ 96  96 135]
    Charge density dimension [ 96  96 135]
```

3. The output above corresponds to the FFT dimensions in our QE calculations. In order to use this charge density and preform charge defect correction we will need to ensure that all our FFT dimensions match these values exactly. Similar to above, before we actually run the charge defect correction we will need to obtain the file `prist.d_tot` (i.e. the total electrostatic potential of the pristine system). Enter the directory `../../Prist/JDFTx` and run the calculation.

```bash
cd ../../Prist/JDFTx
sbatch job
```

Note!!!!! check for the following line in `main.in`

```
    fftbox 96 96 135
```

This line make sures our fft grid our file work in subsequent calculations.

4. Next return to the `../../Sn-Cluster` folder and run the subsequent charged defect calculations in the `Q+*/CCC` folders. Note if we change anything about the FFT grid of the system (i.e. change cutoff, lattice parameter, or directly set a different FFT) then we will need to redo step 3 and recreate a new `prist.d_tot`.

```bash
cd ../../Sn-Cluster
cd Q+1/CCC
mv ../conv_charge.n_* .     # make sure to move the charge density file we created in step 2!
sbatch job
cd ../../Q+2/CCC
mv ../conv_charge.n_* .     # make sure to move the charge density file we created in step 2!
sbatch job
```

Note!!!!! check for the following lines in `main.in`

```
    fftbox 96 96 135
    fix-electron-density conv_charge.$VAR
```

As above we fix the fft grid, but furthermore we will fix electron density (note the absence of `electronic-minimize`!)

5. If you did all the above correctly you should find the following charge cell corrections:

```bash
grep Net Q+*/main.out
```
Produces the output:

```
    Q+1/main.out:	Net correction: 0.00656478 (= EmodelIsolated - EmodelPeriodic + q deltaV)
    Q+2/main.out:	Net correction: 0.01888650 (= EmodelIsolated - EmodelPeriodic + q deltaV)
```


---

Author(s)
------------------------------------
Tyler J. Smart

