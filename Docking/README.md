Description of Directory Contents
=================================

```1wco_N.pdb```- This is the nisin part of the pdb file for nisin downloaded from the proten data bank. It has the molecular structure that corresponds to Nisin A.

```2012_ADTtut.pdf``` - This is the pdf file that explains how to use the Autodock software. Follow the same steps outlined in the pdf but use ```1wco_N.pdb``` in place of ```ind.pdb``` and ```4Y68_A.pdb``` in place of ```hsg1.pdb```.

```4Y68_A.pdb``` - This is one of the 4 NSR monomers downloaded from the Protein Data Bank.

```config.txt``` - This file contains the parameters required by Autodock during the docking procedure for nisin into NSR.

```nisnPROVAL.pdb``` - This is pdb file created from mutating the ```1wco_N.pdn``` file in the Chimera software so that the SER residue and the isoleucine residue change to the proline and valine residues respectively.  


UPDATE: 28/07/18
================

Following the 3 steps in ```2012_ADTtut.pdf``` may result in an error. When performing step 3 and after the HOH* and * have been added a different warning appears saying no atoms have been selected

![screen shot 2018-07-28 at 19 26 33](https://user-images.githubusercontent.com/13021392/43359623-2dc49252-929d-11e8-8a3c-8ff111e5d650.png)
