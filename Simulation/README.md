File Descriptions
=================

This folder contains the files that are used at the start of the MD Simulation. One can generate every other file from these base files by following the steps in the primary README file.

```nisinmol7.pdb``` - This is the starting configuration of the Nisin molecule.  

```4Y68_A.pdb``` - This is the starting configuration of the NSR molecule.  

```DBB.cif``` and  ```DHA.cif``` - These are the ```cif``` files for DBB and DHA which are used to create the frcmod files for DBB and DHA respectively.  

```dbb.mc``` and ```dha.mc``` - These are the ```mc``` file that tell LEaP what atoms to ignore in the DHA and DBB molecules when they part of nisin.  

```min.in``` and ```min2.in``` - These are the input files for use with ```sander``` containing the parameters for energy minimisation.  

```heat.in``` and ```equib*.in``` - These are the input files for use with ```sander``` containing the parameters for the heating and equilbrium stages of the MD simulation.  

```prodser1ns1.in``` - This is the input file for use with ```sander``` containing the parameters for a production stage of 1 nanosecond.

```hbondtraj.ptraj``` - This is the cpptraj file used to generate the raw data for the hbond analysis.  

```trajfiles.atom2.ptraj``` - This is the cpptraj file used to generate the raw data file ```SER236OtoCYS28Cnyl50``` for the distance analysis in Nisin A.

```trajfiles.proval.ptraj``` - This is the cpptraj file used to generate the raw data file ```SER236OtoCYS28Cnylproval50``` for the distance analysis in Nisin PV.

```trajfiles.rmstime50.ptraj``` - This is the cpptraj file used to generate the raw data file ```rmsdnisinovertime``` for the RMSD analysis in Nisin A.  

```trajfiles.rmsdpvovertime.ptraj``` - This the is the cpptraj file used to generate the raw data file ```rmsdpvovertime``` for the RMSD analysis in Nisin PV.  

```trajfiles.rmsd50.ptraj``` - This is the cpptraj file used to generate the raw data file ```perresavg.BB.50.dat``` for the RMSF analysis in Nisin A.

```trajfiles.rmsdpv.ptraj``` - This is the cpptraj file used to generate the raw data file ```perresavg.BB.50proval.dat``` for the RMSF analysis in Nisin PV.   









 
