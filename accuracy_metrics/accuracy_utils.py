# Turn off RDKit logging to ignore the canonicalization errors (https://github.com/rdkit/rdkit/issues/2683)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import pandas as pd
import sys
sys.path.append("../")
# from accuracy_metrics.accuracy_utils import *
# from accuracy_metrics.rxnmapper.rxnmapper import BatchedMapper
from rdkit import Chem

def canonicalize_smiles(smiles: str, remove_atom_number: bool = True):
    # smiles = "".join(smiles.split())
    cano_smiles = None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None
    if mol is None:
        return None
    
    if remove_atom_number:
        for a in mol.GetAtoms():
            if a.HasProp('molAtomMapNumber'):
                a.ClearProp('molAtomMapNumber')
    try:
        cano_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None
        
    return cano_smiles
    
def exact_match(rxn_smiles, prediction_list, n_best):
    accuracy = np.zeros(n_best, dtype=float)

    gt,_,prod = rxn_smiles.strip().split(">")
    gt = canonicalize_smiles(gt)
    prod = canonicalize_smiles(prod)
    
    evaluate_prediction_list = []
    topk_count = 0
    # evaluate predictions
    for prediction in prediction_list[:n_best]:
        prediction = canonicalize_smiles(prediction)
        if prediction is None:
            topk_count += 1
        else:
            # check if unique prediction (issue: https://gitlab.com/mlpds_mit/askcosv2/retro/graph2smiles/-/issues?show=eyJpaWQiOiIyIiwiZnVsbF9wYXRoIjoibWxwZHNfbWl0L2Fza2Nvc3YyL3JldHJvL2dyYXBoMnNtaWxlcyIsImlkIjoxNjg3NTg0OTB9)
            if prediction not in evaluate_prediction_list:
                evaluate_prediction_list.append(prediction)
                if prediction == gt:
                    accuracy[topk_count:] = 1.0
                    break
                topk_count += 1
    return accuracy.tolist()


def superset(rxn_smiles, prediction_list, n_best):
    accuracy = np.zeros(n_best, dtype=float)

    gt,_,prod = rxn_smiles.strip().split(">")
    gt = canonicalize_smiles(gt)
    prod = canonicalize_smiles(prod)
    
    evaluate_prediction_list = []
    topk_count = 0
    gt = set(gt.split("."))
    for prediction in prediction_list[:n_best]:
        prediction = canonicalize_smiles(prediction)
        if prediction is None:
            topk_count += 1
        else:
            prediction = set(prediction.split("."))
            if prediction not in evaluate_prediction_list:
                evaluate_prediction_list.append(prediction)
                if prediction.issuperset(gt):
                    accuracy[topk_count:] = 1.0
                    break
                topk_count += 1
    return accuracy.tolist()


###### Stereochemistry-agnostic ######

def remove_stereochemistry(smiles):
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None
    if mol is None:
        return None
        
    Chem.RemoveStereochemistry(mol)

    try:
        smi_no_stereo = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    except Exception:
        return None
    return smi_no_stereo


def stereochemistry_agnostic(rxn_smiles, prediction_list, n_best):
    accuracy = np.zeros(n_best, dtype=float)
    prod = canonicalize_smiles(rxn_smiles.strip().split(">")[2])
    gt = canonicalize_smiles(remove_stereochemistry(rxn_smiles.strip().split(">")[0]))
    if gt is None or prod is None:
        return accuracy.tolist()

    evaluate_prediction_list = []
    topk_count = 0

    for prediction in prediction_list[:n_best]:
        try:
            prediction = canonicalize_smiles(remove_stereochemistry(prediction))
        except Exception:
            topk_count += 1
            continue
        if prediction is None:
            topk_count += 1
        else:
            if prediction not in evaluate_prediction_list:
                evaluate_prediction_list.append(prediction)
                if prediction == gt:
                    accuracy[topk_count:] = 1.0
                    break
                topk_count += 1
    return accuracy.tolist()

#### Synthon ####

from collections import defaultdict

def remove_unmapped_atoms_and_cap(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # unmapped atom indices (from original mol)
    unmapped = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == 0]
    if not unmapped:
        return smiles

    # count how many bonds each mapped atom loses
    h_to_add = defaultdict(int)

    for atom in mol.GetAtoms():
        if atom.GetIdx() in unmapped:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomMapNum() > 0:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    h_to_add[nbr.GetAtomMapNum()] += int(bond.GetBondTypeAsDouble())

    emol = Chem.RWMol(mol)

    # remove ALL unmapped atoms 
    for idx in sorted(unmapped, reverse=True):
        emol.RemoveAtom(idx)

    # cap mapped atoms w H
    for mapnum, nH in h_to_add.items():
        atom = next(a for a in emol.GetAtoms() if a.GetAtomMapNum() == mapnum)

        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)

        for _ in range(nH):
            h = Chem.Atom(1)
            h_idx = emol.AddAtom(h)
            emol.AddBond(atom.GetIdx(), h_idx, Chem.BondType.SINGLE)

    new_mol = emol.GetMol()

    try:
        Chem.SanitizeMol(
            new_mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                        Chem.SanitizeFlags.SANITIZE_KEKULIZE
            )
    except Exception:
        # this is usually due to the case of metal-complex.
        return None

    return Chem.MolToSmiles(new_mol, canonical=True, isomericSmiles=True)


def synthon(rxn_smiles_rxnmapper, superset_accuracy_list, predicted_mapped_rxn_list, n_best):

    gt_synthon = canonicalize_smiles(remove_unmapped_atoms_and_cap(rxn_smiles_rxnmapper.split(">")[0]))
    prod = canonicalize_smiles(rxn_smiles_rxnmapper.split(">")[2])
    accuracy = np.array(superset_accuracy_list, dtype=float)

    if gt_synthon is None or prod is None:
        return accuracy.tolist()

    gt_synthon_set = set(gt_synthon.split("."))
    topk_count = 0
    evaluate_prediction_list = [] 
    for k in range(n_best):
        if superset_accuracy_list[k] == 1:  # already superset correct
            accuracy[k:] = 1.0
            return accuracy.tolist()
        
        mapped_prediction_rxn = predicted_mapped_rxn_list[k]
        try:
            prediction_synthon = canonicalize_smiles(remove_unmapped_atoms_and_cap(mapped_prediction_rxn.split(">")[0]))
        except Exception:
            topk_count += 1
            continue
        if prediction_synthon is None or prediction_synthon == "":
            topk_count += 1
            continue
        prediction_synthon_set = set(prediction_synthon.split("."))
        if prediction_synthon_set not in evaluate_prediction_list:
            evaluate_prediction_list.append(prediction_synthon_set)
            if prediction_synthon_set == gt_synthon_set:
                accuracy[topk_count:] = 1.0
                break
            topk_count += 1
    return accuracy.tolist()

### Two-Step Superset Accuracy ###

def sort_fragments_by_heavy_atom(fragments):
    fragments = sorted(fragments, key=lambda x: len(Chem.MolFromSmiles(x).GetAtoms()), reverse=True)
    return fragments # first fragment is the biggest fragment which is passed through the trained model again.

def two_step_superset(
    superset_accuracy_list, 
    second_pass_prediction_dict, 
    n_best_two_step=5
):
    """
    Compute the two-step superset accuracy metric.

    This metric evaluates whether, for each of the top-n_best first-pass predictions:
    - The candidate is a superset match to the ground truth precursor (first-pass success).
    - Or, if not, whether any of the top-50 second-pass predictions (combined 
      with small fragments) is a superset match to the ground truth precursor.

    Parameters
    ----------
    superset_accuracy_list : np.ndarray of float
        List of length n_best containing superset accuracy (1/0) for each of the top-n_best first-pass predictions.
    second_pass_prediction_dict : dict
        Dictionary of length n_best. Each value is a dict for a given candidate k containing:
            - 'rxn_smiles': str, the ground-truth reaction SMILES.
            - 'cand_{k+1}_smaller_frags': str or empty, small fragment(s) to append if needed.
            - [f"cand_precursor_{i+1}"]: str, predicted precursor from the second-pass (for i in 0..49).
    n_best_two_step : int (default=5)
        How many top first-pass candidates to consider.

    Returns
    -------
    accuracy : np.ndarray
        Array of floats, shape (n_best,), where each position is 1 if any correct prediction is found at or before this rank.
    """
    if len(second_pass_prediction_dict) != n_best_two_step:
        raise ValueError(f"second_pass_prediction_dict must have {n_best_two_step} predictions for n_best first-pass candidates")
    
    accuracy = np.zeros(n_best_two_step, dtype=float)
    found_correct = False
    for k in range(n_best_two_step):
        if superset_accuracy_list[k] == 1:  # first-pass already superset correct
            accuracy[k:] = 1.0
            return accuracy.tolist()
        else:
            # Evaluate second-pass predictions for candidate k
            smaller_fragment = second_pass_prediction_dict[k][f"cand_{k+1}_smaller_frags"]
            gt = set(canonicalize_smiles(second_pass_prediction_dict[k]["rxn_smiles"].split(">")[0]).split("."))

            prediction_second_pass_list = second_pass_prediction_dict[k][[f"cand_precursor_{i+1}" for i in range(50)]] # consider top-50 predictions from second run.
            for second_pass_pred in prediction_second_pass_list:
                if pd.isna(smaller_fragment) or smaller_fragment == "":
                    prediction = canonicalize_smiles(second_pass_pred)
                else:
                    # only concat if smaller_fragment is a non-empty string
                    if isinstance(smaller_fragment, str) and smaller_fragment.strip() != "":
                        prediction = canonicalize_smiles(str(second_pass_pred) + "." + smaller_fragment)
                    else:
                        prediction = canonicalize_smiles(second_pass_pred)
                if prediction is None:
                    continue
                prediction = set(prediction.split("."))
                if prediction.issuperset(gt):
                    found_correct = True
                    accuracy[k:] = 1.0
                    break
            if found_correct:
                break
    return accuracy.tolist()