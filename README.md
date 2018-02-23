MD-Simulation
=============

This repository holds the codebases for the peptide-protein docking procedure ("Docking" folder) and the codebase for using the Amber suite of prorams ("Simulation" folder) to run a 50 nansecond Molecular Dynamic Simulation. Also here we describe in full the various steps that are required to set up and perform the simulation correctly.  

The molecular modeling simulations consisted of 2 models. The first model used a configuration of residues 22-34 of a nisin mutant where residue 29 is serine and is structurally aligned so that residues 22-34 fit into the tunnel region of the NSR molecule. The second model was also aligned to fit into the tunnel region of the NSR enzyme only this time residue 29 has been mutated to proline.

## Creating the NSR-Nisin pdb file  


To create a pdb file for the NSR-Nisin complex, individual files for both NSR and Nisin are neeed. This requirement can be satified quite easily by using the R package bio3d like so

```{r Initial - version}

> install.packages("bio3d")
> library(bio3d)
> # This code fetches the NSR pdb file and splits it into separate pdb files for each chain > get.pdb("4Y68", URLonly = F, split=TRUE, path="split_chain.test", multi=F)
> # This code fetches the nisin file and spilt the pdb file into individual models
> # with a separate file for each chain
> get.pdb("1wco", URLonly = F, split=TRUE, path="split_chain.nisin", multi=T)
```
For the NSR molecule the pdb file 4Y68 was used and for the nisin molecule the 1wco pdb file was used where Protein Databank Accession numbers for both files came from their respective papers [2,3]

## Docking Procedure
To determine a starting configuration for the NSR-Nisin complex a docking program Autodock for ligand and protein binding was used [4]. This produced 9 possible binding conformations for the NSR-Nisin complex. Out of these 9 states the 7th conformation state was chosen being the state which showed the most favorable interaction between residues 29 of nisin and residues 236-240 of the active site in NSR. In this model the NSR-nisin complex had serine for residue 29. The molecular visualisation package Chimera was then used to mutate residue 29 to proline and residue 30 to valine and the subsequent configuration saved as a different pdb file. Again the docking program Autodock was used to create a starting configuration for the NSR-Nisin complex. In this case the 1st conformation was chosen for the same reasons as previously. Then both chosen conformation states (one with serine, the other with proline) were saved as separate pdb files to be used as starting configurations in the molcular modeling simulation.


## Molecular Modelling 

To run the MD simulation the Amber workflow was used [5]. This consisted of the following steps.

- Prepare the pdb files using LEaP.

- Prepare cif files for nonstandrad residues DBB and DHA.

- Use antechamber forcefield values to fill in for missing parameters.

- Use antechamber to create prmtop and inpcrd files for each pdb file.

- Used sander to run MD simulation from prmtop and inpcrd files.

- Use cpptraj to analyse MD trajectory files for change of distances over time.

- Use cpptraj to determine H-bonds and water mediated interactions.

- Use mmpbsa to calculate binding energies of residues.

In what follows the procedure for the workflow as applied to the NSR-Nisin model with serine is described. The steps for the workflow as applied to the NSR-Nisin model with proline are identical. Here the Amber 16 suite of programs was used [5]. To start the pdb files from the docking procedure step were prepared for use with Amber’s LEaP program [6]. This required running the following Linux command.

```bash

$ pdb4amber -i Nisinmol7.pdb -o gfp.pdb
$ pdb4amber -i 4Y68_A.pdb -o NSR.pdb
$ pdb4amber -i NSR.pdb -o recptor.pdb
$ reduce -Trim recptor.pdb > recptorNoH.pdb

```

Here ```Nisinmol7.pdb``` was the fie created by Autodock for the 7th binding conformation state and ```gfp.pdb``` the output file from the pdb4amber program. And ```4Y68_A.pdb``` was the pdb file for chain A in the NSR protein. The ```pdb4amber``` program changed the residues labled HIS to HIE and indicated that the nonstandard residues dehydroalanine (DHA) and d-alpha-aminobutyric acid (DBB) were not recognised by LEaP. Also the ```reduce``` program with the ```-Trim``` flag strips all hydrogens from the pdb file [7]

To deal with the nonstandard residues dehydroalanine (DHA) and d-alpha-aminobutyric acid (DBB) the respective entries in the RCSB Protein Data Bank were accessed and the respective ```.cif``` files downloaded. The following BASH code was then run using several programs from Amber.

```bash

$ antechamber -fi ccif -i DBB.cif -bk DBB -fo ac -o dbb.ac -c bcc -at amber 
$ prepgen -i dbb.ac -o dbb.prepin -m dbb.mc -rn DBB
$ parmchk2 -i dbb.prepin -f prepi -o dbb.frcmod -a Y -p parm10.dat
$ antechamber -fi ccif -i DHA.cif -bk DHA -fo ac -o dha.ac -c bcc -at amber 
$ prepgen -i dha.ac -o dha.prepin -m dha.mc -rn DHA
$ parmchk2 -i dha.prepin -f prepi -o dha.frcmod -a Y -p parm10.dat
```

The amber program antechamber can read the ```bash .cif``` files of the nonstandard residues and assign partial charges and atom types to the the nonstandard residues based on the bcc charge scheme[8]. This will output ``` .ac``` files which have charge and bonding information for the nonstandard residues. These will then used as input to the prepgen program along with a custom made ``` mc``` file that tells the prepgen program what atoms to ignore from the residue (for peptide bonding). The ``` mc``` file should look like so

```bash
HEAD_NAME N 
TAIL_NAME C 
MAIN_CHAIN CA 
OMIT_NAME OXT 
OMIT_NAME HXT 
PRE_HEAD_TYPE C 
POST_TAIL_TYPE N 
CHARGE 0.0
```

The HEAD NAME and TAIL NAME lines identify the atoms that will connect to the previous and following amino acids, respectively. The MAIN CHAIN lines list the atoms along the chain that connect the head and the tail atoms. The OMIT NAME lines list the atoms in the nonstandard residue that should be removed from the final structure, as they are not present in the intact protein. The PRE HEAD TYPE and POST TAIL TYPE lines let prepgen know what atom types in the surrounding protein will be used for the covalent connection. The CHARGE line gives the total charge on the residue; prepgen will ensure that the charges of the ”omitted” atoms are redistributed among the remaining atoms so that the total charge is correct (i.e., 0 in this case).The prepgen program then outputs ```.prepin``` files which are the inputed into the parmchk2 program that create the ```.frcmod``` files using paramters from the gaff.dat and parm10.dat parameter files. Next the ```.prepin``` and ```.frcmod``` files were read into the LEaP program

```bash

Macintosh-109add6f31eb:111.ROUGH tonyblake$ tleap
-I: Adding /Users/tonyblake/amber16/dat/leap/prep to search path.
-I: Adding /Users/tonyblake/amber16/dat/leap/lib to search path.
-I: Adding /Users/tonyblake/amber16/dat/leap/parm to search path.
-I: Adding /Users/tonyblake/amber16/dat/leap/cmd to search path.

Welcome to LEaP!
(no leaprc in search path)
> source leaprc.protein.ff14SB       #reads in ff14SB forcefield parameters
> set default PBRadii mbondi3        
> loadamberprep dbb.prepin           #reads in prep file for DBB
> loadamberparams dbb.frcmod         #reads in forcefield parameters for DBB
> loadamberprep dha.prepin           #reads in prep file for DHA
> loadamberparams dha.frcmod         #reads in forcefield parameters for DHA
> x=loadPDB gfp.pdb                  #reads in gfp.pdb
> saveamberparm x gfp.prmtop gfp.inpcrd  # creates topology and coordinate files
```

At this point however LEaP began to complain about missing parameters for bond lengths, bond angles and dihedral angles. The different values for these parameters (they are different fior each atom) are used with the amber force field to determine the energies for the NSR-Nisin complex[9].  These values are the values for parameters kb, r0, kθ, θ0, γ, Vn, Aij, Bij which are specified in the ```.frcmod``` and ```.dat``` files. So to overcome the ”missing parameter” issue values from the ```gaff2.dat``` file were used as the values for the missing parameter indicated by LEaP. The underlying cause of the ”missing parameter” issue was LEaP’s inability to recognise the 3 peptide bonds and associated angles and dihedral angles of LYS-DBB, ALA-DBB, and VAL-DHA. Thus the ```extra.frcmod``` file was created to supply these missing values to LEaP. 

```bash

Macintosh-109add6f31eb:003.SANDER tonyblake$ vim extra.frcmod 

Remark line goes here
MASS
NT 14.010        0.530               same as n4
C  12.010        0.616               same as c
O  16.000        0.434               same as o
CT 12.010        0.878               same as c3
CD 12.010        0.360               same as c2
CM 12.010        0.360               same as c2
CX 12.010        0.878               same as c3
N  14.01         0.530

BOND
C -NT  255.5    1.5460             same as c -n4

ANGLE
O-C -NT   69.53     118.830
C -NT-H   44.63     111.120          same as c -n4-hn SOURCE3_SOURCE5
C -NT-CT  62.14     108.760          same as c -n4-c3        SOURCE5
CX-C -NT  64.28     112.260          same as c3-c -n4        SOURCE3
CD-C -N   86.65     124.990           same as c2-c2-nh            SOURCE3
C -NT-CD  63.556    112.580          same as c2-n4-c2

DIHE

CX-C -NT-H    1    1.025       180.000          -2.
O-C -NT-CT    4   10.000       180.000           2.
CX-C -NT-CT   9    1.400         0.000           3.
O-C -NT-H     1    2.500       180.000          -2.
C -NT-CD-CM   6    0.000       180.000           3.
CX-C -NT-CD   9    1.400         0.000           3.
O-C -NT-CD    4   10.000       180.000           2.
C -NT-CD-C    6    0.000        80.000           3.
CM-CD-C -N    4   26.600       180.000           2.
NT-CD-C -N    4   26.600       180.000           2.

IMPROPER

NONBON
```

There was one other small issue which seems to be a bug with the LEaP program. In some instances after the extra.frcmod was read into leap and the command run for creating the topology and coordinate files LEaP would still not recognise the parameters for 1 or 2 dihedral angles

```bash

> loadamberparams extra.frcmod         
Loading parameters: ./extra.frcmod
Reading force field modification type file (frcmod)
Reading title:
Remark line goes here
> saveamberparm x gfp.prmtop gfp.prmcrd
Checking Unit.
WARNING: The unperturbed charge of the unit: 2.000000 is not zero.

 -- ignoring the warning.

Building topology.
Building atom parameters.
Building bond parameters.
Building angle parameters.
Building proper torsion parameters.
 ** No torsion terms for  O-C-NT-CD
 ** No torsion terms for  CX-C-NT-CD
Building improper torsion parameters.
old PREP-specified impropers:
 <DBB 2>:  CA   +M   C    O   
 <DBB 4>:  CA   +M   C    O   
 <DHA 12>:  C    CB   CA   N   
 <DHA 12>:  CA   HB1  CB   HB2 
 <DHA 12>:  CA   +M   C    O   
 total 30 improper torsions applied
 5 improper torsions in old prep form
Building H-Bond parameters.
Incorporating Non-Bonded adjustments.
Parameter file was not saved.
```
To solve this issue one only needs to move the the lines in the extra.frcmod corresponding to those 2 dihedral angles with missing parameters to the top of the dihedral angle parameter list. So if the above example had to be changed then the corrected version would look like so.

```bash

DIHE

O-C -NT-CD    4   10.000       180.000           2.
CX-C -NT-CD   9    1.400         0.000           3.
C -NT-CD-C    6    0.000        80.000           3.
CX-C -NT-H    1    1.025       180.000          -2.
O-C -NT-CT    4   10.000       180.000           2.
CX-C -NT-CT   9    1.400         0.000           3.
O-C -NT-H     1    2.500       180.000          -2.
C -NT-CD-CM   6    0.000       180.000           3.
CM-CD-C -N    4   26.600       180.000           2.
NT-CD-C -N    4   26.600       180.000           2.

IMPROPER

NONBON
```

Then LEaP is able to read all of the nisin molecule (```gfp.pdb```) and will create the corresponding topology (```gfp.prmtop```) and coordinate files (```gfp.prmcrd```). When this happens LEaP will output the following onto the screen

```bash
> loadamberparams extra.frcmod         
Loading parameters: ./extra.frcmod
Reading force field modification type file (frcmod)
Reading title:
Remark line goes here
> saveamberparm x gfp.prmtop gfp.prmcrd
Checking Unit.
WARNING: The unperturbed charge of the unit: 2.000000 is not zero.

 -- ignoring the warning.

Building topology.
Building atom parameters.
Building bond parameters.
Building angle parameters.
Building proper torsion parameters.
Building improper torsion parameters.
old PREP-specified impropers:
 <DBB 2>:  CA   +M   C    O   
 <DBB 4>:  CA   +M   C    O   
 <DHA 12>:  C    CB   CA   N   
 <DHA 12>:  CA   HB1  CB   HB2 
 <DHA 12>:  CA   +M   C    O   
 total 30 improper torsions applied
 5 improper torsions in old prep form
Building H-Bond parameters.
Incorporating Non-Bonded adjustments.
Not Marking per-residue atom chain types.
Marking per-residue atom chain types.
  (Residues lacking connect0/connect1 - 
   these dont have chain types marked:

	res	total affected

	CLYS	1
	NLYS	1
  )
 (no restraints)
```

Then to carry out the rest of the LEaP procedure for the NSR enzyme, solvating nisin and NSR, combining nisin and NSR, and adding counterions the following commands were issued.

```bash
> source leaprc.water.tip3p
> solvateBox x TIP3PBOX 10.0 
> saveamberparm x gfp.wat.prmtop gfp.wat.prmcrd 
> REC=loadPDB recptor.pdb
> NSR=loadPDB recptorNoH.pdb
> saveamberparm NSR nsr.prmtop nsr.prmcrd
> solvateBox NSR TIP3PBOX 10.0
> saveamberparm NSR nsr.wat.prmtop nsr.wat.prmcrd
> LIG=loadPDB gfp.pdb
> saveamberparm LIG test.prmtop test.prmcrd
> PROT=loadPDB recptorNoH.pdb
> saveamberparm PROT test2.prot.prmtop test2.prot.prmcrd
> COM = combine {PROT LIG}
> saveamberparm COM com.prmtop com.prmcrd
> solvateBox COM TIP3PBOX 10.0
> saveamberparm COM com.wat.prmtop com.wat.prmcrd
> addIons x Cl- 0  
> saveamberparm x gfp.neutral.prmtop gfp.neutral.prmcrd
> addIons LIG Cl- 0
> saveamberparm LIG gfp.neutral2.prmtop gfp.neutral2.prmcrd
> saveamberparm x gfp.wat.neutral.prmtop gfp.wat.neutral.prmcrd
> addIons PROT Cl- 0 
> saveamberparm PROT nsr.neutral2.prmtop nsr.neutral2.prmcrd 
> addIons COM Cl- 0 
> saveamberparm COM com.neutral2.prmtop com.neutral2.prmcrd
> saveamberparm COM com.wat.neutral2.prmtop com.wat.neutral2.prmcrd
> addIons NSR Cl- 0
> saveamberparm NSR nsr.wat.neutral2.prmtop nsr.wat.neutral2.prmcrd 
> LIGNEUT=loadPDB gfp.pdb
> saveamberparm LIGNEUT test3.prmtop test3.prmcrd
> PROTNEUT=loadPDB recptorNoH.pdb 
> COMNEUT=combine {PROTNEUT LIGNEUT}  
> addIons COMNEUT Cl- 0
> save COMNEUT com.neutral0.prmtop com.neutral0.prmcrd
> saveamberparm COMNEUT com.neutral0.prmtop com.neutral0.prmcrd
```
These commands produced many files that needed to be used with ```sander```, ```cpptraj``` and ```MMPBSA.py```. Next the sander program was used to carry out the molecular simulation. To use this program an input file needed to be created that would inform the program of the physical processes to simu- late. Noting the good results of a previous study ?, the same parameters were employed. Firstly the system was minimisied in a two stage process. The input file for the first minimisation stage looks like this.

