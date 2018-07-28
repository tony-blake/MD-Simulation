Description of Directory Contents
=================================

```1wco_N.pdb```- This is the nisin part of the pdb file for nisin downloaded from the proten data bank. It has the molecular structure that corresponds to Nisin A.

```2012_ADTtut.pdf``` - This is the pdf file that explains how to use the Autodock software. Follow the same steps outlined in the pdf but use ```1wco_N.pdb``` in place of ```ind.pdb``` and ```4Y68_A.pdb``` in place of ```hsg1.pdb```.

```4Y68_A.pdb``` - This is one of the 4 NSR monomers downloaded from the Protein Data Bank.

```config.txt``` - This file contains the parameters required by Autodock during the docking procedure for nisin into NSR.

```nisnPROVAL.pdb``` - This is pdb file created from mutating the ```1wco_N.pdn``` file in the Chimera software so that the SER residue and the isoleucine residue change to the proline and valine residues respectively.  


UPDATE: 28/07/18
================

Following the steps in Exercies 1 in ```2012_ADTtut.pdf``` may result in an error. When performing step 3 and after the HOH* and * have been added a different warning appears saying no atoms have been selected

![screen shot 2018-07-28 at 19 26 33](https://user-images.githubusercontent.com/13021392/43359623-2dc49252-929d-11e8-8a3c-8ff111e5d650.png)

If this happens close Autodock (not saving the project) and open up Chimera. You can use Chimera to perform Autodock Prep (essentially Exercise 1 as would have been performed in the auto dock procedure). After starting Chimera open the protein file (4Y68_A.pdb). From the pull down menu go to "tools" in the menu bar and select "Surface Binding Analysis" and then "Dock-Prep". Then click "OK" for the next 4 pop-up windows. Then save the file as ```4Y68_A.mol2```. (For a visual guide to these steps see the ```034-chimera_vina.pdf``` document). Then close Chimera. Then open up the ```4Y68_A.mol2``` file in Chimera. Then save it as ```protein.pdb``` Then open up Autodock and proceed from Exercise 2 of ```2012_ADTtut.pdf```. In Exercise 3 select the new protein file you created from Chimera Dock-prep procedure (```protein.pdb```) as the molecule to be selected and save it as a pdbqt file as per the instructions in the exercise. Then follow the rest of the exercises to complete the docking procedure.
