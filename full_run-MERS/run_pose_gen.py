from dask.distributed import Client
from rdkit import Chem
import traceback
import dask
import tqdm
from dask.distributed import LocalCluster
import fegrow
import time
import pathlib

from rdkit.Chem import rdFMCS

def find_best_n_matches(target_ligand, ref_ligands, n_top):
    matched_mols = []
    for ref_lig in ref_ligands:
        mcs = rdFMCS.FindMCS(
            [target_ligand, ref_lig],
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            atomCompare=rdFMCS.AtomCompare.CompareAnyHeavyAtom,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            maximizeBonds=False,
            timeout=1,
        )
        matched_mols.append((mcs.numAtoms, mcs.numBonds, ref_lig, mcs.smartsString))
    matched_mols.sort(key=lambda x: (x[0], x[1]), reverse=True)
    return matched_mols[:n_top]

def transfer_coordinates(reference_ligand: Chem.Mol, template_ligand: Chem.Mol) -> Chem.Mol:
    """
    Transfer the coordinates from the reference to the template ligand.

    Args:
        reference_ligand: The ligand we want to generate the conformers for.
        template_ligand: The ligand whose coordinates should be used as a reference.

    Returns:
        The template ligand with atom positions set to the reference for overlapping atoms.

    """
    matches = reference_ligand.GetSubstructMatch(template_ligand)
    if not matches:
        raise RuntimeError(
            f"A core fragment could not be extracted from the reference ligand using core smarts {Chem.MolToSmarts(template_ligand)}"
        )

    ref_conformer: Chem.Conformer = reference_ligand.GetConformer(0)
    template_conformer = Chem.Conformer()
    for i, atom_match in enumerate(matches):
        ref_atom_position = ref_conformer.GetAtomPosition(atom_match)
        template_conformer.SetAtomPosition(i, ref_atom_position)
    template_ligand.AddConformer(template_conformer, assignId=True)
    return template_ligand

@dask.delayed
def pose_ligand(target_ligand, reference_ligands, protein_file, name) -> dict:
    # returns the posed ligand and the ligand used as the template
    # for the input ligand generate a pose using the best mcs match as the refernce ligand
    best_matches = find_best_n_matches(target_ligand=target_ligand, ref_ligands=reference_ligands, n_top=1)
    for close_match in best_matches:
        # try each match till we can generate a pose
        # now create the template molecule with the correct coords
        mcs_mol = Chem.MolFromSmarts(close_match[-1])
        mcs_mol = transfer_coordinates(reference_ligand=close_match[2], template_ligand=mcs_mol)
        # make the fegrow mol and use this as the template
        fe_mol = fegrow.RMol(target_ligand)
        fe_mol._save_template(mcs_mol)
        # generate coords we found 300 to work well in asap-alchemy
        try:
            fe_mol.generate_conformers(num_conf=300)
        except ValueError:
            print('Failed to embed trying next match')
            continue
        # if we remove all of the conformers keep cycling through the list to see if this helps
        fe_mol.remove_clashing_confs(protein=protein_file)
        # optimise with ani2x
        fe_mol.optimise_in_receptor(receptor_file=protein_file, ligand_force_field="openff")
        # sort the conformers so the lowest energy is at id 0
        fe_mol.sort_conformers()
        if fe_mol.GetNumConformers() != 0:
            break
        else:
            print("Removed all conformers trying next match")
    # now purge all but the lowest in energy
    fe_mol = fegrow.RMol(fe_mol, confId=0)
    return {'posed_mol':fe_mol, 'template_mol':fegrow.RMol(close_match[2]), 'name': name}



def main():
    workers = 20
    client = Client(LocalCluster(threads_per_worker=1, n_workers=workers))
    # create the cluster
    print(f"Client created {client}")

    # load all of the training ligands
    training_ligands = []
    lig_file = pathlib.Path('new_core.sdf')
    rdkit_mol = Chem.MolFromMolFile(lig_file.as_posix(), removeHs=True)
    training_ligands.append(rdkit_mol)

    # load the test molecules
    test_mols = []
    supplier = Chem.SmilesMolSupplier("smiles-test-MERS.txt")
    for lig in supplier:
        lig = Chem.AddHs(lig)
        test_mols.append(lig)

    ref_ligands = dask.delayed(training_ligands)

    # build a list of tasks to submit
    for_submission = [
        pose_ligand(
            target_ligand=test_lig,
            reference_ligands=ref_ligands,
            protein_file=str(pathlib.Path('rec_final.pdb').absolute()),
            name=f'ligand_{i}'
        )
        for i, test_lig in enumerate(test_mols)
    ]
    submitted = client.compute(for_submission)

    output_path = pathlib.Path('outputs')
    output_path.mkdir(exist_ok=True)
    with tqdm.tqdm(total=len(submitted), desc="Scoring molecules...", ncols=80) as pbar:
        while len(submitted) > 0:
                to_remove = []
                for job in submitted:
                    if not job.done():
                        continue

                    # remove the job
                    to_remove.append(job)
                    pbar.update(1)


                    try:
                        mol_data = job.result()
                        name = mol_data['name']
                        out_dir = output_path.joinpath(name)
                        out_dir.mkdir(exist_ok=True, parents=True)
                        posed_ligand = mol_data['posed_mol']
                        temp_lig = mol_data['template_mol']
                        posed_ligand.to_file(out_dir.joinpath("best_pose.sdf").as_posix())
                        temp_lig.to_file(out_dir.joinpath("template_lig.sdf").as_posix())
                        
                    except Exception as E:
                        traceback.print_exc()

                # remove jobs before next round
                for job in to_remove:
                    submitted.remove(job)

                time.sleep(5)

    print("All molecules scored")


if __name__ == '__main__':
    main()