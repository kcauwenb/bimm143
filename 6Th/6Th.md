6Th
================

Section 1: In silico docking of drugs to HIV-1 protease
-------------------------------------------------------

We want to produce a protein-only PDB file and drug-only PDB file.

``` r
library(bio3d)
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

Q1: What is the name of the two non protein resid values in this structure? What does resid correspond to and how would you get a listing of all reside values in this structure?

HOH and MK1

``` r
prot <- trim.pdb(hiv, "protein")
protein <- atom.select(hiv, "protein",value=T)
lig <- trim.pdb(hiv, "ligand")
ligand <- atom.select(hiv, "ligand",value=T)
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

Q2: Can you locate the binding site visually? Note that crystal structures normally lack hydrogen atoms, why?

No.

Section 2: Docking ligands into HIV-1 protease
----------------------------------------------

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

Section 3: Exploring the conformational dynamics of proteins
------------------------------------------------------------
