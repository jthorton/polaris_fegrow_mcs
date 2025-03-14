# Polaris pose challenge

Submission to the polaris binding pose challenge.

## Polaris FEGrow-MCS

[FEGrow](https://github.com/cole-group/FEgrow) is an open-source software package which aims to automate the building of a 
congeneric series of ligands within a target protein building on experimental structural information. 

Typically, this involves the user selecting a single core molecule who's structure is preserved within in a series as the reference and new r-groups and linkers are built 
and refined while holding the core fixed. This however can result in the core being very small to ensure it overlaps with an entire series of diverse ligands. In this case 
it might be better to have the core structure automatically derived for each target ligand to be built as the maximum common substructure (mcs) between the reference and target to
ensure we use as much of the experimental structural information as possible. 

We can then extend this idea to cases where we have multiple reference molecules, now rather than having to extract a core for each reference we can use the MCS overlap to 
automatically extract a core ligand. 

## Protocol

After the core structure was found an ensemble of up to 300 conformations is generated using the ETKDG algorithm in RDKit, with restraints on the core structure. Conformers with 
heavy clashes with the protein are then removed and the remaining ensemble is optimised using OpenMM while holding the protein positions fixed. AMBER FF14SB is used to 
parameterise the protein and OpenFF-Sage with am1-bcc charges from ambertools is used to represent the non-bonded terms of the target ligand, while the ANI-2x ML potential is used 
to treat the intramolecular energy of the ligand. After energy minimising each of the conformers with the hybrid ML/MM method the lowest energy conformation is extracted as the best pose.

### SARS

For each of the test ligands the 10 largest, sorted by mcs bonds and atoms, overlapping mcs cores were extracted as core templates to grow the ligand. The ligand was grown from each core 
in order of decreasing mcs overlap until a successful set of optimised poses could be built. The lowest energy pose was then submitted. 

The structures were then manually inspected and a list of 14 ligands were highlighted as having collapsed structures where the ML potential had failed, these we re-run using OpenFF-Sage for the intramolecular energies 
of the ligand as well. 

Some molecules still failed to build and had to be replaced with ligands from a default FEGrow run see [here](https://github.com/cole-group/polaris-fegrow).

### MERS

As we only have a single training ligand for this set it was used to extract an MCS core for all test set ligands. The core was first rebuilt using FEGrow due to a strange bend 
in the vector from the pyridine core which prevented all structures from being built. Some structures still consistently failed to be built this way and so they 
were replaced with ligands from a default FEGrow run see [here](https://github.com/cole-group/polaris-fegrow).