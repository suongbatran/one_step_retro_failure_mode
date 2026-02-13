import sys
sys.path.append("../")


import json
import logging
from typing import Set
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

RS_TAG = ["R", "S", "None"]
RS_TAG_DICT = {rs: i for i, rs in enumerate(RS_TAG)}


def analyze_molecule(smiles: str) -> Set[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # Return empty sets and 0 for invalid SMILES
        return set(), set(), 0
        
    atom_types_mapped = set()
    atom_types_unmapped = set()
    num_unmapped_atoms = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum():
            atom_types_mapped.add(atom.GetAtomicNum())
        else:
            atom_types_unmapped.add(atom.GetAtomicNum())
            num_unmapped_atoms += 1
    return atom_types_mapped, atom_types_unmapped, num_unmapped_atoms

unspecified_atom = Chem.ChiralType.CHI_UNSPECIFIED
unspecified_bond = Chem.BondStereo.STEREONONE

def has_atom_stereochemistry(atom):
    return atom.GetChiralTag() != unspecified_atom

def has_bond_stereochemistry(bond):
    return bond.GetStereo() != unspecified_bond

def has_atom_stereo_from_smiles(smiles):

    return "@" in smiles

def has_bond_stereo_from_smiles(smiles):
    return "/" in smiles or "\\" in smiles

def get_rs_tag(atom):

    return atom.GetPropsAsDict().get("_CIPCode", "None")

sp3 = Chem.rdchem.HybridizationType.SP3
double = Chem.rdchem.BondType.DOUBLE

def check_stereochemistry_changes(smiles):

    # Note: Checks for stereochemistry changes, keeps a record of the hybridization/bond types as well

    try:
        reactants, _, products = smiles.split(">")
        r_mol = Chem.MolFromSmiles(reactants)
        p_mol = Chem.MolFromSmiles(products)

        has_p_atom_stereo = any(
            has_atom_stereochemistry(atom) 
            for atom in p_mol.GetAtoms()
        )
        has_p_bond_stereo = any(
            has_bond_stereochemistry(bond)
            for bond in p_mol.GetBonds()
        )

        atom_changes = {}
        bond_changes = {}

        # Map atom map numbers to atoms in reactants
        r_atoms = {atom.GetAtomMapNum(): atom for atom in r_mol.GetAtoms() if atom.GetAtomMapNum() > 0}
        p_atoms = {atom.GetAtomMapNum(): atom for atom in p_mol.GetAtoms() if atom.GetAtomMapNum() > 0}

        # Check for atom chirality changes
        for map_num, p_atom in p_atoms.items():
            r_atom = r_atoms.get(map_num)
            if not r_atom:
                continue

            p_rs_tag = get_rs_tag(p_atom)
            r_rs_tag = get_rs_tag(r_atom)

            if (
                p_rs_tag != r_rs_tag
                # and r_atom.GetHybridization() == sp3
                # and p_atom.GetHybridization() == sp3
            ):
            # if p_atom.GetChiralTag() != r_atom.GetChiralTag() and r_atom.GetHybridization() == sp3:
                atom_changes[map_num] = (
                    (r_rs_tag, r_atom.GetHybridization() == sp3),
                    (p_rs_tag, p_atom.GetHybridization() == sp3)
                )

        # Create a map from (sorted atom map numbers) to bond in reactants
        def get_bond_key(bond, mol):
            begin_map = mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomMapNum()
            end_map = mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomMapNum()
            if begin_map == 0 or end_map == 0:
                return None
            return "-".join(map(str, sorted([begin_map, end_map])))

        r_bonds = {
            get_bond_key(bond, r_mol): bond
            for bond in r_mol.GetBonds()
            if get_bond_key(bond, r_mol)
        }

        # Compare bond stereo
        for p_bond in p_mol.GetBonds():
            bond_key = get_bond_key(p_bond, p_mol)
            if not bond_key or bond_key not in r_bonds:
                continue
            r_bond = r_bonds[bond_key]
            if (
                p_bond.GetStereo() != r_bond.GetStereo() 
                # and r_bond.GetBondType() == double
                # and p_bond.GetBondType() == double
            ):
                bond_changes[bond_key] = (
                    (r_bond.GetStereo(), r_bond.GetBondType() == double),
                    (p_bond.GetStereo(), p_bond.GetBondType() == double)
                )
        atom_changes = json.dumps(atom_changes) if atom_changes else ""
        bond_changes = json.dumps(bond_changes) if bond_changes else ""

        return has_p_atom_stereo, has_p_bond_stereo, atom_changes, bond_changes
    except Exception as e:
        

        # Default to checking the smiles if error in creating mol obkect
        r, _, p = smiles.split(">")
        has_p_atom_stereo = has_atom_stereo_from_smiles(p)
        has_p_bond_stereo = has_bond_stereo_from_smiles(p)

        logging.error(f"Error checking stereochemistry inversion: {e} with smiles: {smiles}, "
            f"product check result - atom: {has_p_atom_stereo}, bond: {has_p_bond_stereo}")

        return has_p_atom_stereo, has_p_bond_stereo, "", ""
