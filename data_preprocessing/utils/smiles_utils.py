from rdkit import Chem
from rdkit.Chem import rdChemReactions

def remove_parentheses(reaction_smiles):
    if reaction_smiles:
        try:
            rxnmol = rdChemReactions.ReactionFromSmarts(reaction_smiles)
            rmols = rxnmol.GetReactants()
            amols = rxnmol.GetAgents()
            pmols = rxnmol.GetProducts()

            reactants = [Chem.MolToSmiles(rmol) for rmol in rmols]
            reagents = [Chem.MolToSmiles(amol) for amol in amols] 
            products = [Chem.MolToSmiles(pmol) for pmol in pmols]

            reaction_smiles = '.'.join(reactants) + '>' + '.'.join(reagents) + '>' + '.'.join(products)
            return reaction_smiles
        except Exception as e:
            print(f"Error processing SMILES: {reaction_smiles}\n{e}")
            return None 
    return None 

def has_any_mapping(mol):
    return any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()) if mol else False

def has_double_mapped_atoms(mol_list):
    mapping_counts = {}
    for mol in mol_list:
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                mapping_counts[map_num] = mapping_counts.get(map_num, 0) + 1
    return any(count > 1 for count in mapping_counts.values())

def canonicalize_smiles(smiles):
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        can_smiles = Chem.MolToSmiles(mol, canonical=True)
        if ('->' in can_smiles) or ('<-' in can_smiles):
            can_smiles = can_smiles.replace('->', '')
            can_smiles = can_smiles.replace('<-', '')
        return can_smiles
    else:
        return None

        
def num_unmapped_atoms(mol):
    count_unmapped = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            count_unmapped +=1
    return count_unmapped

def has_at_least_1C(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return any(atom.GetSymbol() == "C" for atom in mol.GetAtoms())

def remove_mapping_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        nomap_smiles = Chem.MolToSmiles(mol)
        return canonicalize_smiles(nomap_smiles)

def remove_mapping_rxnsmiles(rxn_smiles):
    reactants, _, product  = rxn_smiles.split('>')
    reactants_smi = reactants.split('.')
    products_smi = product.split('.')
    
    new_reactants_smi = []
    for smi in reactants_smi:
        new_smi = remove_mapping_smiles(smi)
        if new_smi:
            new_reactants_smi.append(new_smi)
    
    new_products_smi = []
    for smi in products_smi:
        new_smi = remove_mapping_smiles(smi)
        if new_smi:
            new_products_smi.append(new_smi)

    if len(new_reactants_smi) != 0 and len(new_products_smi) != 0:
        return '.'.join(sorted(new_reactants_smi))+ '>>' + '.'.join(sorted(new_products_smi))

def get_neighbor_info(smiles):
    """Extracts a mapping of atom index → sorted list of its neighbor mapping numbers."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return {
            atom.GetAtomMapNum(): sorted([nbr.GetAtomMapNum() for nbr in atom.GetNeighbors() if nbr.GetAtomMapNum() != 0])
            for atom in mol.GetAtoms() if atom.GetAtomMapNum() != 0
        }
    else:
        return None
