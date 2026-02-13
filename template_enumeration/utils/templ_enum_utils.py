import os
import logging
from rdkit import Chem, RDLogger
from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRun
from utils.chem_utils import canonicalize_smiles

RDLogger.DisableLog('rdApp.*')

def has_passed_filter(prec, template):
    smiles_list = prec.split('.')
    if template.get("intra_only") and len(smiles_list) > 1:
        # Disallowed intermolecular reaction
        return False
    if template.get("dimer_only") and (
            len(set(smiles_list)) != 1 or len(smiles_list) != 2
    ):
        # Not a dimer
        return False
    return True


def apply_template_to_product(template, product, prod=None, candidate_small_frag=None):
    """
    Applies the template to the product and returns a list of canonicalized precursor SMILES.
    If prod (rdchiralReactants object) is not provided, it will be created from product.
    """
    try:
        smarts = template['reaction_smarts']
        rxn = rdchiralReaction(smarts)
        if prod is None:
            prod = rdchiralReactants(product)
        precs = rdchiralRun(rxn, prod, combine_enantiomers=False)

        precs = [
            prec for prec in precs 
            if has_passed_filter(prec, template)
        ]

        if candidate_small_frag:
            return [canonicalize_smiles(f"{prec}.{candidate_small_frag}") for prec in precs]
        else:
            return [canonicalize_smiles(prec) for prec in precs]
    except Exception as e:
        logging.error(f"Error applying template {template['reaction_smarts']}: {e}")
        return []