MD-Simulation
=============

This repository holds the codebases for the peptide-protein docking procedure ("Docking" folder) and the codebase for using the Amber suite of prorams ("Simulation" folder) to run a 50 nansecond Molecular Dynamic Simulation. Also here we describe in full the various steps that are required to set up and perform the simulation correctly.  

The molecular modeling simulations consisted of 2 models. The first model used a configuration of residues 22-34 of a nisin mutant where residue 29 is serine and is structurally aligned so that residues 22-34 fit into the tunnel region of the NSR molecule. The second model was also aligned to fit into the tunnel region of the NSR enzyme only this time residue 29 has been mutated to proline.

### Creating the NSR-Nisin pdb file  


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


### Molecular Modelling 

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

At this point however LEaP began to complain about missing parameters for bond lengths, bond angles and dihedral angles. The different values for these parameters (they are different fior each atom) are used with the amber force field to determine the energies for the NSR-Nisin complex[9]. 

