# How to reproduce real time DMD pipeline from scratch?
## JDFTx
[JDFTx](https://jdftx.org) is the first principles code that is used to do SCF, phonon, and wannier calculations before running 
the DMD calculations. The JDFTx source code is available from [github](https://github.com/shankar1729/jdftx/tags). You can download 
it to directly to the cluster, for example, the following command,
   
```wget https://github.com/shankar1729/jdftx/archive/refs/tags/v1.7.0.tar.gz```
  
will download the version 1.7 of the code. The tar ball can be unarchived and decompressed using the following command, 

```tar -xvzpf v1.7.0.tar.gz```. 

Once uziped change to the main directory, ```/home/PATH_TO_CODE/jdftx-1.7.0```, replace PATH_TO_CODE with actual path relative to your home directory.
From inside the main directory if you issue command ```ls```, you will see something very similar to the following:   
```Dockerfile  fluid1D  generate_tarball.sh  jdftx```.   

Here you can issue the command, ```mkdir build```, to create the build directory 
where the compiled code or executables would be saved. Now using command ```cd build``` change to the build directory. Inside the build directory
using a text editor of your choice create a text file (```make.sh``` for example, name is really not important) and save the following script in it and issue ```sh make.sh```. 

```cmake
# for UCSC lux cluster using gcc compilers

module load openmpi/gcc/64/1.10.7

#Configure:
CC=mpicc CXX=mpicxx cmake \
   -D EnableMKL=yes \
   -D MKL_PATH="/cm/shared/apps/intel/mkl" \
   -D EnableScaLAPACK=yes \
   -D GSL_PATH="/cm/shared/apps/gsl/2.6" \
   -D ForceFFTW=yes \
   -D FFTW3_PATH="/cm/shared/apps/fftw/fftw-3.3.8" \
   -D EnableLibXC=yes \
   -D LIBXC_PATH="/home/YOUR_USERNAME/Programs/lib/libxc-6.0.0" \
   -D EnableProfiling=yes \
   ../jdftx/

make -j8
```

```cmake
# for stampede2 using intel compilers

module load intel/18.0.2 impi/18.0.2 gsl/2.6

#Configure:
CC=mpicc CXX=mpicxx cmake \
   -D EnableMKL=yes \
   -D MKL_PATH="$TACC_MKL_DIR" \
   -D EnableScaLAPACK=yes \
   -D GSL_PATH="$TACC_GSL_DIR" \
   -D ForceFFTW=yes \
   -D FFTW3_PATH="$TACC_FFTW3_PATH" \
   -D EnableLibXC=yes \
   -D LIBXC_PATH="/home/YOUR_USERNAME/Programs/lib/libxc-6.0.0" \
   -D EnableProfiling=yes \
   ../jdftx/

make -j8
```

There are several bits in the above script which are specific to the HPC platform that you are working with. For example
```module load openmpi/gcc/64/1.10.7``` tells cmake which compilers to use. In this case it instructors to use ```gcc```. Moreover the exact module 
```openmpi/gcc/64/1.10.7``` may not be available. You may want to use ```module avail``` command to check which ```openmpi``` or ```impi``` modules
are available and change accordingly. 

The other platform dependent bits are values of ```MKL_PATH```, ```GSL_PATH```, and ```FFTW3_PATH``` and would need to be changed.
So far the four applications, ```gcc``` or ```icc```, ```MKL```, ```GSL```, and ```FFTW3```, we have discussed are very common to scientific computing
and most HPC platforms will have them already. You would not need to install them. Just find out how to correctly access them. The ```LibXC``` library, 
however, is not usually available on most platforms and you may need to build it yourself and then provide its path to the cmake script above. 

Although a basic compilation can be done without some of the optional libraries but the DMD calculations are computationally expensive and would 
need an optimally customized build. The above script does just that. Also these are condensed instructions to get one started quickly. For comprehensive
information on how to build a fully customized JDFTx code please consult the [JDFTx](https://jdftx.org/Compiling.html) documentation from the developers.

**Note:** For DMD calculations a separate code FynWann is also needed. FynWann uses JDFTx as a library and would give a compile time error if original JDFTx is not
compiled with the profiling option enabled, ```EnableProfiling=yes```.

### Building Libxc
Libxc is a library of exchange-correlation functionals for density functional theory. You can get the library from [Gitlab](https://gitlab.com/libxc/libxc/-/releases), for example, using the ```wget``` command from within the cluster;   
```wget https://gitlab.com/libxc/libxc/-/archive/6.0.0/libxc-6.0.0.tar.gz```.   
Next decompress it with,    
```tar -xvfz libxc-6.0.0.tar.gz```.
After that change to the main directory ```/PATH_TO_LIBRARY/libxc-6.0.0```. If you do not see the ```configure``` file in the main directory, you can create one by issuing:   

```autoreconf -i```.

After that you can issue the following commands in order to build the directroy:

```
CFLAGS="$CFLAGS -std=gnu99" ./configure --enable-shared --prefix="/PATH_TO_LIBRARY/libxc-6.0.0"
make
make check
make install
```   
**Note 1:** The libxc has some code which requires to be compiled with C99 standard, that is why without the compiler flag ```-std=gnu99``` libxc would not compile.   
**Note 2:** JDFTx is built both as an executable and a library (shared object). That is why it is required that libxc is built as a shared object as well otherwise you will have problem compiling JDFTx with libxc. That is why the configure option ```--enable-shared``` is very important.  

### FFTW3

```./configure --enable-shared --enable-mpi --enable-threads --prefix=/software/groups/ping_group/shared/libs/fftw-3.3.10/build```

### FeynWann

FeynWann is another code which works with JDFTx (and needs compiled JDFTx to compile) to perform important energy calculations and inialiazations for
the DMD calculations. The code is not publically available yet. If you are working on the project you will have private access to it. Download the code and save it at the same location where your JDFTx ```build``` direcotry is. Make a new directory, say, ```build-FeynWann```, and change to it. Inside this directroy save the following script to a file ```make-FeynWann.sh```. 

```cmake
#!/bin/bash
# for building FeynWann on lux
module load openmpi/gcc/64/1.10.7

CC=mpicc CXX=mpicxx cmake \
 -D JDFTX_BUILD="../build" \
 -D JDFTX_SRC="../jdftx" \
 -D ForceFFTW=yes \
 -D FFTW3_PATH="/PATH_TO_LIBRARY/fftw-3.3.8" \
 -D EnableMKL=yes \
 -D MKL_PATH="/PATH_TO_MKL/mkl" \
 -D EnableScaLAPACK=yes \
 -D GSL_PATH="/PATH_TO_LIBRARY/gsl/" \
 -D EnablePETSc=yes \
 -D PETSC_PATH="/PATH_TO_LIBRARY/petsc-3.18.4-build" \
 -D MPISafeWrite=no \
../FeynWann

make -j4
```
### PETSc 
**Download a release tatball instead of git clone as it may lead to problems with pre-generated Fortran-stubs.**
The PETSc library optionally used by FeynWann may not be available on some systems. Here is how to install. Download it from [here](https://petsc.org/release/install/download/#recommended-obtain-release-version-with-git). Inside the main directory make a new directory called ```build``` and run the configure command with following options,

```./configure --prefix=/PATH_TO_LIBRARY/petsc-3.18.4/build --with-blaslapack-dir=/PATH_TO_MKL/mkl --with-mpi-dir=/PATH_TO_OPENMPI/openmpi/gcc/64/1.10.7```.

The configure command prints the next command to run, if there is no problem with the configure options. Run those commands in sequence one after the other to complete the PETSc installation. 

### GSL
GNU Scientific Library provides a wide range of mathematical routines. Here is how to install it. 
Inside the main directory make a new diroctory called ```build``` and run the following command,

```./configure --prefix=/software/groups/ping_group/shared/libs/gsl-2.7.1/build```

Then issue ```make``` and ```make install``` in sequence. 

The above setup is required for first of the two parts of DMD simualtions. This part helps obtain the electronic and phonon structure at the Kohn-Sham level and necessary initializations for the DMD calculations. At this stage it is best to run through an example calculation to obtain the electronic and phonon structure using JDFTx. 
## GaAs
### Electronic band structure
We take the example of GaAs for this purpose and first compute its electronic structure. There are two algorithms to find the converged electronic ground state. We first perform SCF and then use the variational minimize - this helps obtain a fully coverged ground state relatively quickly. For the SCF calculations you need input file - very later convenience - split into two files, ```common.in``` and ```scf.in```. 

```
# save this in common.in

lattice face-centered Cubic 10.6829
ion-species Ga_nv3_nocorecorr.upf
ion-species As_nv5_nocorecorr.upf
elec-cutoff 17

ion Ga 0.00 0.00 0.00  0
ion As 0.25 0.25 0.25  0

elec-n-bands 34
converge-empty-states yes
spintype spin-orbit

elec-ex-corr mgga-x-scan mgga-c-scan
```

```
# save this in scf.in

include common.in

kpoint-folding 24 24 24    #Use a Brillouin zone mesh

initial-state totalE.$VAR

electronic-scf energyDiffThreshold 1e-8

dump-name totalE.$VAR
dump Init Symmetries
dump End State
```
These two ```.in``` files specify the structure and other simulations parameters. Along with these two files, you also need to specify either the path to the directory which has pseudopotentials or save the pseudopotential files in the same directory. The pseudpotentials used in this example are available from here, ![Ga](pseudos/Ga_nv3_nocorecorr.upf) and ![As](pseudos/As_nv5_nocorecorr.upf). Here is typical script to run the above SCF calculation with JDFTx. 
```
#!/bin/bash
#SBATCH -p cpuq
#SBATCH --account=cpuq
#SBATCH -N 8
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=8
#SBATCH -J jdftx

module load openmpi/gcc/64/1.10.7

pwd; hostname; date

echo "Running program on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS total tasks"
echo  "with each node getting $SLURM_NTASKS_PER_NODE tasks."

MPICMD="mpirun -np $SLURM_NTASKS"
DIRJ="/PATH_TO_EXECUTABLE/JDFTx/jdftx-1.7.0/build"
${MPICMD} ${DIRJ}/jdftx -i scf.in > scf.out
```
After doing the SCF convergence, we repeat the ground state covergence with a different, slower but more reliable, algorithm called ```electronic minimize``` to do the final convergence. It uses the SCF converged state and further converges it. Here is the input for this step,

```
# save the following in totalE.in
include common.in

kpoint-folding 24 24 24    #Use a Brillouin zone mesh

initial-state totalE.$VAR

electronic-minimize energyDiffThreshold 1e-11

dump-name totalE.$VAR
dump End State EigStats BandEigs ElecDensity Vscloc
```
And repeat the calculations by just changing this line in the above batch script,

```${MPICMD} ${DIRJ}/jdftx -i totalE.in > totalE.out```

Once this calculation has converged, the electronic band structure can be calculated. Before doing the 
electronic band structure calculations we must list the high symmetry points in the Brillouin zone to decide
the path along which we want to calculate the band structure. Save the following high symmetry points for GaAs in a file, for
example, ```bandstruct.kpoints.in```, 

```
kpoint 0.000 0.000 0.000     Gamma
kpoint 0.000 0.500 0.500     X
kpoint 0.250 0.750 0.500     W
kpoint 0.500 0.500 0.500     L
kpoint 0.000 0.000 0.000     Gamma
kpoint 0.375 0.750 0.375     K
```
Now to generate k-points along this path run the following, 

```/PATH_TO_EXECUTABLE/JDFTx/jdftx-1.7.0/jdftx/scripts/bandstructKpoints bandstruct.kpoints.in 0.01 bandstruct```

This script will generate two files, ```bandstruct.kpoints``` and ```bandstruct.plot```. The first file has the k-points 
aloing the specified path and the second file has the script to plot the calculated bands. 
Now you can use the following input to run jdftx for bandstructure calculations,

```
# Save the following in bandstruct.in
include common.in
include bandstruct.kpoints

dump-name bandstruct.$VAR
dump End BandEigs Spin

electronic-minimize energyDiffThreshold 1e-11
fix-electron-potential totalE.$VAR
```
Once this is done, you can open the ```bandstruct.plot``` file to modify it a little to make a nice looking band structure plot, 
here is the example, 

```gnuplot
#!/usr/bin/gnuplot -persist
set term pngcairo
set output "GaAs-band.png"

set xtics ( "Gamma" 0,  "X" 142,  "W" 213,  "L" 284,  "Gamma" 458,  "K" 642 )
unset key
nRows = real(system("awk '$1==\"kpoint\" {nRows++} END {print nRows}' bandstruct.kpoints"))
nCols = real(system("wc -c < bandstruct.eigenvals")) / (8*nRows)
formatString = system(sprintf("echo '' | awk 'END { str=\"\"; for(i=0; i<%d; i++) str = str \"%%\" \"lf\"; print str}'", nCols))
set xzeroaxis               #Add dotted line at zero energy
set ylabel "E - VBM [eV]"   #Add y-axis label
set yrange [*:10]           #Truncate bands very far from VBM

plot for [i=1:nCols] "bandstruct.eigenvals" binary format=formatString u 0:((column(i) - 0.180239)*27.21) w l lw 2
```
The above gnuplot scritp will make a figure that looks like the one shown below, 
![Bands](figs/GaAs-band-4x4x4.png)


### Phonon Dispersion
Similarly the phonon dispersion for the same system can be calculated. Although one needs a bigger simulation cell for 
a proper phonon calculation. A supercell made of 4x4x4 unit cells is reasonable cell for this system. Use the following 
jdftx input to calculate force matrix. Remember to change the executable in the batch script: ```${MPICMD} ${DIRJ}/phonon -i phonon.in > phonon.out```.

```
include totalE.in          #Full specification of the unit cell calculation
initial-state totalE.$VAR  #Start from converged unit cell state
dump-only                  #Don't reconverge unit cell state

phonon supercell 4 4 4     #Calculate force matrix in a 4x4x4 supercell
```

After the above calculation, the phonon dispersion can be caculated using the following python scrip, 

```python
# save the following to PhononDispersion.py:
import numpy as np
from scipy.interpolate import interp1d

#Read the phonon cell map and force matrix:
cellMap = np.loadtxt("phonon.phononCellMap")[:,0:3].astype(np.int)
forceMatrix = np.fromfile("phonon.phononOmegaSq", dtype=np.float64)
nCells = cellMap.shape[0]
nModes = int(np.sqrt(forceMatrix.shape[0] / nCells))
#Read the k-point path:
kpointsIn = np.loadtxt('bandstruct.kpoints', skiprows=2, usecols=(1,2,3))
nKin = kpointsIn.shape[0]
#--- Interpolate to a 10x finer k-point path:
nInterp = 10
xIn = np.arange(nKin)
x = (1./nInterp)*np.arange(1+nInterp*(nKin-1)) #same range with 10x density
kpoints = interp1d(xIn, kpointsIn, axis=0)(x)
nK = kpoints.shape[0]

#Calculate dispersion from force matrix:
#--- Fourier transform from real to k space:
forceMatrixTilde = np.tensordot(np.exp((2j*np.pi)*np.dot(kpoints,cellMap.T)), forceMatrix, axes=1)
#--- Diagonalize:
omegaSq, normalModes = np.linalg.eigh(forceMatrixTilde)

#Plot phonon dispersion:
import matplotlib.pyplot as plt
meV = 1e-3/27.2114
plt.plot(np.sqrt(omegaSq)/meV)
plt.xlim([0,nK-1])
plt.ylim([0,None])
plt.ylabel("Phonon energy [meV]")
#--- If available, extract k-point labels from bandstruct.plot:
try:
    import subprocess as sp
    kpathLabels = sp.check_output(['awk', '/set xtics/ {print}', 'bandstruct.plot']).split()
    kpathLabelText = [ label.split('"')[1] for label in kpathLabels[3:-2:2] ]
    kpathLabelPos = [ nInterp*int(pos.split(',')[0]) for pos in kpathLabels[4:-1:2] ]
    plt.xticks(kpathLabelPos, kpathLabelText)
except:
    print ('Warning: could not extract labels from bandstruct.plot')
#plt.xticks ( "Gamma" 0,  "X" 71,  "W" 107,  "L" 143,  "Gamma" 230,  "K" 322 )
plt.xticks ([nInterp*0,nInterp*142,nInterp*213,nInterp*284,nInterp*458, nInterp*642],[r"$\Gamma$", "X", "W", "L", r"$\Gamma$", "K"])
plt.savefig("phononDispersion-4x4x4.png", format="png", bbox_inches="tight")
plt.show()
```
The above script will produce a graph that looks like the one below. 

![Phonon](figs/phononDispersion.png)

### Wannier Band Structure

In the next stage we calculate maximally localized wannier functions (MLWFs). The wannier calculation uses two files (quantities) 
```totalE.Zeff``` containing **Born effective charges** and ```totalE.epsInf``` containing **optical dielectric tensor** which need to be computed externally, for example, using Quantum Espresso code. Another important piece of input for wannier calculation is the initial guess for the wannier centers. There is no unique way to do obtain this guess. One way of doing it is using randomly generated wannier centers. Here is a python script which randomly generates wannier function centers and saves them to the file ```rand_wann-centers.dat```,

```python
#!/usr/bin/env python
import os
import random

random.seed()

f = open("rand_wann-centers.dat","w")

for i in range(8):
  x=random.random()-0.5
  y=random.random()-0.5
  z=random.random()-0.5
  f.write("wannier-center Gaussian %10.6f %10.6f %10.6f 1.7 %s \n" % (x,y,z,"sUp"))
  f.write("wannier-center Gaussian %10.6f %10.6f %10.6f 1.7 %s \n" % (x,y,z,"sDn"))
````

Once these quantities have been computed externally and saved in the same directory, the following 
jdftx input can be used for wannier calculations,

```
# save the following to wannier.in
include totalE.in

wannier \
        loadRotations yes \
        innerWindow  -0.11 0.318 \
        outerWindow  -0.11 0.7 \
        saveMomenta yes \
        saveSpin yes \
        phononSupercell 4 4 4 \
        polar yes

wannier-initial-state totalE.$VAR
wannier-dump-name wannier.$VAR

wannier-minimize \
    energyDiffThreshold 1e-11 \
    nIterations 50000

wannier-center Gaussian  -0.299074  -0.449527  -0.280677 1.7 sUp
wannier-center Gaussian  -0.299074  -0.449527  -0.280677 1.7 sDn
wannier-center Gaussian   0.314828  -0.125568   0.279016 1.7 sUp
wannier-center Gaussian   0.314828  -0.125568   0.279016 1.7 sDn
wannier-center Gaussian  -0.239254   0.440158  -0.419614 1.7 sUp
wannier-center Gaussian  -0.239254   0.440158  -0.419614 1.7 sDn
wannier-center Gaussian  -0.414207   0.484450  -0.019686 1.7 sUp
wannier-center Gaussian  -0.414207   0.484450  -0.019686 1.7 sDn
wannier-center Gaussian  -0.159616   0.379131   0.208221 1.7 sUp
wannier-center Gaussian  -0.159616   0.379131   0.208221 1.7 sDn
wannier-center Gaussian   0.042114   0.268505  -0.023349 1.7 sUp
wannier-center Gaussian   0.042114   0.268505  -0.023349 1.7 sDn
wannier-center Gaussian  -0.363296  -0.067796  -0.166150 1.7 sUp
wannier-center Gaussian  -0.363296  -0.067796  -0.166150 1.7 sDn
wannier-center Gaussian   0.000263   0.058292   0.442630 1.7 sUp
wannier-center Gaussian   0.000263   0.058292   0.442630 1.7 sDn
```
The initial guess for wannier functions generated using python script has been appended to the above input file. 
Now do the wannier calculation by changing the executable in the batch script:

```${MPICMD} ${DIRJ}/wannier -i wannier.in > wannier.out```.

After that run the following python script to obtain the wannier band structure, 

```python
#Save the following to WannierBandstruct.py:
import numpy as np
from scipy.interpolate import interp1d

#Read the MLWF cell map, weights and Hamiltonian:
cellMap = np.loadtxt("wannier.mlwfCellMap")[:,0:3].astype(int)
Wwannier = np.fromfile("wannier.mlwfCellWeights")
nCells = cellMap.shape[0]
nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
Wwannier = Wwannier.reshape((nCells,nBands,nBands)).swapaxes(1,2)

#--- Get k-point folding from totalE.out:
for line in open('totalE.out'):
    if line.startswith('kpoint-folding'):
        kfold = np.array([int(tok) for tok in line.split()[1:4]])
kfoldProd = np.prod(kfold)
kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])

#print('nCells = ', nCells, "Wwannier.shape = ", Wwannier.shape[0],"nBands = ", nBands, "nK = ", kfoldProd)
#--- Read reduced Wannier Hamiltonian and expand it:
Hreduced = np.fromfile("wannier.mlwfH", dtype=complex).reshape((kfoldProd,nBands,nBands)).swapaxes(1,2)
iReduced = np.dot(np.mod(cellMap, kfold[None,:]), kStride)
Hwannier = Wwannier * Hreduced[iReduced]

#Read the band structure k-points:
kpointsIn = np.loadtxt('bandstruct.kpoints', skiprows=2, usecols=(1,2,3))
nKin = kpointsIn.shape[0]
#--- Interpolate to a 10x finer k-point path:
xIn = np.arange(nKin)
x = 0.1*np.arange(1+10*(nKin-1)) #same range with 10x density
kpoints = interp1d(xIn, kpointsIn, axis=0)(x)
nK = kpoints.shape[0]

#Calculate band structure from MLWF Hamiltonian:
#--- Fourier transform from MLWF to k space:
Hk = np.tensordot(np.exp((2j*np.pi)*np.dot(kpoints,cellMap.T)), Hwannier, axes=1)
#--- Diagonalize:
Ek,Vk = np.linalg.eigh(Hk)
#--- Save:
np.savetxt("wannier.eigenvals", Ek)
```
And finally use the following gnuplot script to plot the wannier bands and the DFT bands together for comparison, 

```gnuplot
#!/usr/bin/gnuplot -persist
set term pngcairo
set output "GaAs-DFTvsWannier.png"

set xtics ( "Gamma" 0,  "X" 142,  "W" 213,  "L" 284,  "Gamma" 458,  "K" 642 )
unset key
nRows = real(system("awk '$1==\"kpoint\" {nRows++} END {print nRows}' bandstruct.kpoints"))
nCols = real(system("wc -c < bandstruct.eigenvals")) / (8*nRows)
formatString = system(sprintf("echo '' | awk 'END { str=\"\"; for(i=0; i<%d; i++) str = str \"%%\" \"lf\"; print str}'", nCols))
set xzeroaxis               #Add dotted line at zero energy
set ylabel "E - VBM [eV]"   #Add y-axis label
set yrange [*:10]           #Truncate bands very far from VBM

VBM = +0.180239  #HOMO from totalE.eigStats
eV = 1/27.2114   #in Hartrees
plot for [i=1:nCols] "bandstruct.eigenvals" binary format=formatString u 0:((column(i) - VBM)/eV) every 14 w p lc rgb "black", \
     for [i=1:10] "wannier.eigenvals" u (0.1*$0):((column(i)-VBM)/eV) w l lw 2
```
The above script will produce a figure similar to the following, 

![Wannier](figs/GaAs-DFTvsWannier.png)
    
This completes the first stage of the calculations. The next stage is the real time density matrix dynamics (DMD). There are two parts to this last stage, first 
is the Lindblad inialization and the second one is DMD run. **It import to note that at the moment the lindblad inialization code does not run on 
lux cluster.**

Using the followig input for Lindblad inialization with all the earlier results saved in the sub-directory ```Wannier``` saved relative to this input file. 

```
scissor         0.424240

NkMult          8
#NkzMult        1
dmuMin          1.43
dmuMax          1.43
Tmax            300
pumpOmegaMax    1.47
pumpTau         50
probeOmegaMax   2.0

ePhMode         DiagK
ePhOnlyElec     1
ePhDelta        0.005
nEphDelta       1000
writeU          1
```

Then change the executable in the batch script to run this input. Here is a batch script to run this on kairay cluster, 

```
#!/bin/bash
##SBATCH -p debug
#SBATCH -N 16
#SBATCH -t 99:00:00
#SBATCH --ntasks-per-node=16
#SBATCH -J lindbladInit


module use /home/jxu153/modulefiles
module load myopenmpi-4.0.2_gcc-4.8.5

MPICMD="mpirun -np $SLURM_NTASKS"
DIRJ="/export/data/share/jxu/jdftx_codes/jdftx-202209/build"
DIRF="/export/data/share/jxu/jdftx_codes/jdftx-202209/build-FeynWann"

${MPICMD} ${DIRF}/lindbladInit_for-DMD-4.5.6/init_for-DMD -i lindbladInit.in > lindbladInit.out
```
The inialization code will save all its output in a sub-directory named, ```ldbd_data```. 
And now one can run the final step of the calculation, which is the DMD. The DMD input is, 

Save the following in the file ```param.in```. 

```
#DEBUG = 1
restart = 0
#compute_tau_only = 1
code = jdftx
alg_linearize = 1
print_along_kpath = 1
kpath_start1 = 0,0,0
kpath_end1 = 0.375,0.75,0.375
alg_only_eimp = 0
alg_eph_sepr_eh = 1
alg_eph_need_elec = 1
alg_eph_need_hole = 0
alg_scatt = lindblad
alg_set_scv_zero = 1
alg_Pin_is_sparse = 0
alg_sparseP = 0
#alg_ddmdteq = 1

#pumpMode = coherent
#pumpA0 = 0.002
pumpPoltype = LC
pumpE = 1.44
pumpTau = 100
probePoltype1 = LC
probePoltype2 = RC
probeEmin = 0.8
probeEmax = 2.0
probeDE = 0.005
probeTau = 100

Bzpert = 0.1

t0 = 0
tend = 1e5
tstep = 100
tstep_pump = 10
freq_measure_ene = 10
de_measure = 5e-4
degauss_measure = 2e-3

#alg_ode_method = euler
#ode_hstart = 0.01
#ode_hmin = 0.01
#ode_epsabs = 1e-04

mu = 1.4
carrier_density = 2e16

#alg_phenom_relax = 1
tau_phenom = 0.11
bStart_tau = 0
bEnd_tau = 4

scrMode = medium
srcFormula = RPA
dynamic_screening = static
epsilon_background = 12.9

impurity_density = 2e16
impMode = model_ionized
#partial_ionized = 1
E_impurity = 1.43
Z_impurity = 1
g_impurity = 2
freq_update_eimp_model = 10

eeMode = Pee_update
#eeMode = Pee_fixed_at_eq
freq_update_ee_model = 10
```
The script to run this is, 

```
#!/bin/bash
##SBATCH -p debug
#SBATCH -N 4
#SBATCH -t 99:00:00
#SBATCH --ntasks-per-node=16
#SBATCH -J dm

MPICMD="mpirun -np $SLURM_NTASKS"
DIRDM="/export/data/share/jxu/denmat-codes/bin"

$MPICMD $DIRDM/denmat_dynm_v4.5.6 > out
```
This calculations will produce several output files for analysis. Here is a python script which
does a curve fitting to the spin evolution curve to get spin-relaxtion time, 

```python
#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import scipy.optimize

fin="sz_elec_tot.out"
tstart = 0 # ps
tend = np.inf # ps. if < 0, means +inf
Bfield = 0 # in-plane (for Sz) magnetic field in Tesla

t_data = np.loadtxt(fin,usecols=(0,))
nt = t_data.shape[0]
t_data = t_data * 2.4188843265857e-5 # / 40000. # a.u. to ps
mask = (t_data >= tstart) & (t_data <= tend)
t_data = t_data[mask]
print("Fitting range: [",t_data[0],", ",t_data[-1],"] ps")
t_data = t_data - t_data[0]
s_data = np.loadtxt(fin,usecols=(1,))[0:nt][mask]

# ==== theoretical parameter values ====
beta = 1 / 100. # initial guess of rate in 1/ps
Bfield = Bfield / 2.3505175675871e5 # convert to a.u.
omega = Bfield / 2.4188843265857e-5
phi = 0 # phase
params = beta, omega, phi

# ==== model ====
def decay(t, beta):
  return s
def residuals1(args, t, s):
  return s - decay(t, *args)
def residuals2(args, t, s):
  return s - decay_cosine(t, *args)

# ==== fitting using curve_fit ====
if omega == 0:
  params_cf, _ = scipy.optimize.curve_fit(decay, t_data, s_data)
  params_lsq, _ = scipy.optimize.leastsq(residuals1, beta, args=(t_data, s_data))
else:
  params_cf, _ = scipy.optimize.curve_fit(decay_cosine, t_data, s_data)
  params_lsq, _ = scipy.optimize.leastsq(residuals2, params, args=(t_data, s_data))

print("Global exponential fit by two ways:")
if omega == 0:
  print("cf: tau = ",1/params_cf[0]," ps")
  print("lsq: tau = ",1/params_lsq[0]," ps")
else:

if omega == 0:
  s_fit = decay(t_data, params_lsq[0])
else:
  s_fit = decay_cosine(t_data, params_lsq[0], params_lsq[1], params_lsq[2])
np.savetxt("sfit.out", np.transpose([t_data, s_data, s_fit]))

#log fit
s_data = np.log(np.abs(s_data))
nt = t_data.shape[0]
t1 = np.zeros(nt-1)
for it in range(2,nt+1):
  fit = np.polyfit(t_data[it-2:it],s_data[it-2:it],1)
  t1[it-2] = -1. / fit[0] # ps

print("Log fit of time-resolved spin lifetime:")
if nt-1 > 20:
  print(t1[0:10])
  print(t1[nt-11:nt-1])
else:
  print(t1)
np.savetxt("tau_t.dat",np.transpose([t_data[0:nt-1],t1]))
```

And finally here a plot of the time resolved spin lifetime, 

![Spin Lifetime](figs/spinlifetime.png)





