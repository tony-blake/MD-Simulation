R code for producing plots in figure X
======================================

In this directory are the R files for taking the raw data generated from Amber during the molecular dynamics.

```All.UU.avg.all10proval.dat``` - Raw data file produced from the ```hbonds.ptraj``` file for Nisin PV. It shows values for the hbond occupancies for all residues in the NSR-Nisin complex. This file was then converted a different format and renamed ```hbondprovalsnisinPROVAL2``` for input into R.  

```Bindingenergy.R``` - This is the R code to create the plot showing hbond occupances for specific residues.  

```FINAL_DECOMP_MMPBSA50PROVAL.xls``` - Excel file showing the binding energies for Nisin-PV calculated from the ```MMPBSA.py``` program.  
```FINAL_DECOMP_MMPBSAser50.xlsx``` - Excel file showing the binding energies for Nisin-PV calculated from the ```MMPBSA.py``` program.  

```SER236OtoCYS28Cnyl50``` - Raw data file generated from the ```trajfiles.atom.ptraj```. It shows the value of the distance from the between the cleavage point nitrogen in Residue 28 in Nisin A and the Sidechain O in SER236 in NSR over 50 nanoeconds.  

```SER236OtoCYS28Cnylproval50``` - Raw data file generated from the ```trajfiles.atom2.ptraj```. It shows the value of the distance from the between the cleavage point nitrogen in Residue 28 in Nisin PV and the Sidechain O in SER236 in NSR over 50 nanoeconds.  

```distance.R``` - This is the R code to create the plot showing  value of the distance from the between the cleavage point nitrogen in Residue 28 in Nisin A and Nisin PV and the Sidechain O in SER236 in NSR over 50 nanoeconds.  

```hbonds.R``` - This is the R code to create the plot showing all hbond occupancies > 20% in Nisin.

```rmsd.R``` - This is the R code to create the plot showing the bar chart of Rmsf values averaged over time for each of the residues in Nisin. It also include the code to plot the RMSD over the duration of the simulation.

```rmsdnisinovertime``` - This is the raw data generated from the ```trajfiles.rmstime50.ptraj``` file. It shows values for the RMSD at each time point for Nisin A.

```rmsdpvovertime``` - This is the raw data generated from the ```trajfiles.rmsdpv.ptraj``` file. The orgiginal output file was ```rms_vs_time.BB.50proval.dat```. To be inputed inot R the file had to be converted to another format and renamed ```rmsdpvovertime```. It shows values for the RMSD at each time point for Nisin PV.` 





