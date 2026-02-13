import numpy as np
import pandas as pd
from typing import List, Tuple, Optional, Set, Dict, Any, Union

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator


RDLogger.DisableLog('rdApp.*')

def canonicalize_smiles(smiles: str, remove_atom_number: bool = True) -> str:
    """Adapted from Molecular Transformer"""
    smiles = "".join(smiles.split())
    cano_smiles = ""

    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        if remove_atom_number:
            [a.ClearProp('molAtomMapNumber') for a in mol.GetAtoms()]

        cano_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        # Sometimes stereochem takes another canonicalization... (just in case)
        mol = Chem.MolFromSmiles(cano_smiles)
        if mol is not None:
            cano_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

            # Additional processing added (modified from the template-relevance repo)
            if ('->' in cano_smiles) or ('<-' in cano_smiles):
                cano_smiles = cano_smiles.replace('->', '')
                cano_smiles = cano_smiles.replace('<-', '')
                
    return cano_smiles
    
def canonicalize_smarts(smarts: str) -> str:
    templ = Chem.MolFromSmarts(smarts)
    if templ is None:
        logging.info(f'Could not parse {smarts}')
        return smarts

    canon_smarts = Chem.MolToSmarts(templ)
    if '[[se]]' in canon_smarts:            # strange parse error
        canon_smarts = smarts

    return canon_smarts

