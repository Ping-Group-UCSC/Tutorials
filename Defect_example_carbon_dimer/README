How to calculate the properties of defect (e.g. CBCN)?

----------------------------------------------------------
Build a Defect Structure
----------------------------------------------------------
1. Build a supercell
 (a) The supercell we use in this tutorial is 6x6 hBN. Then replace one B and one adjacent N atoms with two C atoms.
     1) use VESTA to open hBN unit cell file and expand the unit cell to 6x6x1
     2) replace one B and one adjacent N atoms with two C atoms
     3) save as p1 format

2. Relax the geometry
 (a) Enter folder cbcn
 (b) $ sbatch job

----------------------------------------------------------
Electronic Structure
----------------------------------------------------------
1. Calculate the electrostatic potential in the vaccum for aligning the electronic structure with respect to vaccum
 (a) $ cd job_vacuum
 (b) $ sbatch job_ppavex
 (c) $ vi avg.out
 (d) The first column is z positions. 
     The second column is the electrostatic potential at each z position. 
     The largest value (V_vac) of the electrostatic potential is the vacuum that we want to align the electronic band structure.

2. Calculate VBM (valence band maximum) and CBM (conduction band minimum) from pristine hBN
 (a) Enter folder pristine_hBN
 (b) $ sbatch job
 (c) $ vi relax.out
 (d) $ grep highest relax.out

3. Plot the electronic structure
 (a) Read the single-particle eigenvalues of CBCN in the spinup and spindown channel.
 (b) Plot the eigenvalues together with VBM and CBM of pristine hBN with respect to the vacuum electrostatic potential
     by substracting V_vac from the eigenvalues, VBM and CBM.

4. Plot wavefunctions
 (a) $ cd job_ppx
 (b) $ sbatch job_ppx
 (c) Folder wfc is where the wavefunctions are stored

5. Analyze symmetry of the wavefunctions and label the state according to the symmetry group which CBCN belongs to
 (a) use VESTA to view one wavefunction
 (b) figure out the irreducible representation that the wavefunction transforms like
 (c) label the wavefunction-corresponding state by the lowercase of the irreducible representation


-----------------------------------------------------------
Charged Defect Formation Energy
-----------------------------------------------------------
coming soon ...


reference:
    https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.110.095505
    https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.102.016402
    https://pubs.rsc.org/en/content/articlepdf/2015/nr/c4nr06376g
    https://journals.aps.org/prb/pdf/10.1103/PhysRevB.86.045112
    https://pubs.aip.org/aip/jcp/article/146/10/104109/195119/First-principles-electrostatic-potentials-for

-----------------------------------------------------------
Nonradiative Recombination and Lifetime
-----------------------------------------------------------
1. Optimize ground state geometry (relax-gs)

2. Optimize excited state geometry (relax-cdftup1) with constrained-DFT

3. ZPL = E(cdftup1) - E(gs)

4. Do linear extrapolation for the geometry between the ground state and excited state

5. Run the linear extrapolation files to obtain the potential energy surfaces of the ground state and excited state

6. Run nonrad code to calculate electron-phonon coupling with nonradiative lifetime.


reference:
    https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.081407
    https://journals.aps.org/prb/pdf/10.1103/PhysRevB.90.075202
    https://pubs.aip.org/aip/jap/article/119/18/181101/142538/Tutorial-Defects-in-semiconductors-Combining


-----------------------------------------------------------
Optical Properties and Radiative Lifetime
-----------------------------------------------------------
1. Optimize defect supercell with nk 3x3x1, which has been tested and is close to converge for GW/BSE

2. Run nscf to expand the number of bands (nbnd) to 1800, which is sufficient to converge GW/BSE@PBE

3. Check the pseudopotential files are in the QE wavefunction output folder

4. Run p2y to convert wavefunctions from QE format to Yambo format

5. Calculate G0W0 (one-shot GW) to obtain quasiparticle correction (includes electron correlation) for single-particle levels

6. Calculate eps0 to obtain the inverse of the dielectric screening matrix that will be input for screened exchange potential

7. Use G0W0 quasiparticle energy and eps0 to calculate BSE, that is absorption spectrum 


reference:
    https://pubs.rsc.org/en/content/articlelanding/2019/TC/C9TC02214G#!divAbstract
    https://pubs.aip.org/aip/apl/article/115/21/212101/37562/Carbon-dimer-defect-as-a-source-of-the-4-1-eV


-----------------------------------------------------------
Photoluminescence
-----------------------------------------------------------
1. Optimize ground state geometry (relax-gs)

2. Optimize excited state geometry (relax-cdftup1) with constrained-DFT

3. Run phonon calculation with the ground state geometry

4. Run PL code


reference:
    https://iopscience.iop.org/article/10.1088/1367-2630/16/7/073026/pdf
    https://journals.aps.org/prmaterials/supplemental/10.1103/PhysRevMaterials.6.L042201/SI.pdf
    https://pubs.acs.org/doi/pdf/10.1021/acs.jpca.0c07339
