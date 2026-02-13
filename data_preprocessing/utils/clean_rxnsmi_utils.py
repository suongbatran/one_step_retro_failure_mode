# Turn off RDKit logging to ignore the canonicalization errors (https://github.com/rdkit/rdkit/issues/2683)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit import Chem
from smiles_utils import *

def clean_rxn_smiles(rxn_smiles, min_num_product_heavy_atoms = 5):
    
    """"
    Takes in an atom-mapped reaction SMILES string and returns a new SMILES of the reaction 
    with matched atom mapping style (i.e. only atoms in both reactants and products are mapped)
    
    With the following cleaning steps: 
        1) Reagents moved to reactants then
            a) chemicals in both reactants and products moved to reagents
            b) remove product(s) if the number of heavy atoms in product is less than min_num_product_heavy_atoms
        
        2) Atom mapping checking: 
            Return None if 
            a) reagents or products have 2 atoms having same mapping number ()
            b) the same atom mapping atom number in the reactants and products are the different atoms
            c) Remove atom mapping for atoms that are not in both reactants and products then:
            - remove unmapped reactants
            - remove unmapped products
        
        3) Retain reaction with single product only
    """
    
    if rxn_smiles:
        # split rxn_smiles
        parts = rxn_smiles.split('>')
        if len(parts) == 3: 
            reactants, reagents, products = parts
        else:
            return None
    else:
        return None
    
    

    # 1) move reagents to reactants (there are cases when reagents are also mapped)
    reactants = reactants + '.' + reagents if reagents else reactants
    # try canonicalize 
    reactants_smi = [canonicalize_smiles(smi) for smi in reactants.split('.')]
    products_smi = [canonicalize_smiles(smi) for smi in products.split('.')]
    
    # a) move items in both sides to reagents
    common = set(reactants_smi) & set(products_smi)
    reactants_smi = [canonicalize_smiles(smi) for smi in reactants_smi if smi not in common]
    products_smi = [canonicalize_smiles(smi) for smi in products_smi if smi not in common]
    
    reagents = [reagent for reagent in common]

    # make sure reactants_smi and products_smi are not empty
    if len(reactants_smi) == 0 or len(products_smi) == 0:
        return None
    if None in reactants_smi or None in products_smi:
        return None      
    
    # b) remove product(s) if the number of heavy atoms in product is less than min_num_product_heavy_atoms

    product_mols = []
    for product in products_smi:
        product_mol = Chem.MolFromSmiles(product)
        if product_mol.GetNumHeavyAtoms() >= min_num_product_heavy_atoms:
            product_mols.append(product_mol)
    
    reactant_mols = []
    for reactant in reactants_smi:
        reactant_mol = Chem.MolFromSmiles(reactant)
        reactant_mols.append(reactant_mol)

    # 2) Atom mapping checking
    # a) check if double-mapped  
    reactant_atom_map = {}
    reactant_mapping = []
    if not has_double_mapped_atoms(reactant_mols):
        for mol in reactant_mols:
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() != 0:
                    reactant_mapping.append(atom.GetAtomMapNum())
                    reactant_atom_map[atom.GetAtomMapNum()] = atom.GetAtomicNum()
    else:
        return None
                
    product_atom_map = {}
    product_mapping = []
    if not has_double_mapped_atoms(product_mols): 
        for mol in product_mols:
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() != 0:
                    product_mapping.append(atom.GetAtomMapNum())
                    product_atom_map[atom.GetAtomMapNum()] = atom.GetAtomicNum()
    else:
        return None
    
    # b) check if same atom identity got mapped
    # get atom mapping that are in both products and reactants
    matched_mapping_set = set(product_mapping) & set(reactant_mapping)

    for atom_map in matched_mapping_set:
        if reactant_atom_map[atom_map] != product_atom_map[atom_map]:
            return None
    
    # c) remove atom mapping for atoms that are not in both reactants and products
    # then remove atom mapping if the atom mapping is not in the matched set
    new_product_mols = []
    for mol in product_mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() not in matched_mapping_set:
                atom.SetAtomMapNum(0)
        # a) only take mapped products
        if has_any_mapping(mol):
            new_product_mols.append(mol)
    
    new_reactant_mols =[]
    
    for mol in reactant_mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() not in matched_mapping_set:
                atom.SetAtomMapNum(0)
        
        # b) move unmapped reactants to reagents
        if has_any_mapping(mol):
            new_reactant_mols.append(mol)
        else:
            reagents.append(canonicalize_smiles((Chem.MolToSmiles(mol))))

    # 3) Retain reaction with signle product only
    if len(new_product_mols) != 1:
        return None
        
    # turn to strings list
    new_products_smi = canonicalize_smiles(Chem.MolToSmiles(new_product_mols[0]))
    
    new_reactants_smi=[]
    for mol in new_reactant_mols:
        new_reactants_smi.append(canonicalize_smiles(Chem.MolToSmiles(mol)))

    if (new_products_smi is not None) and all(smi is not None for smi in new_reactants_smi):
    
        # return clean rxn_smiles with reagents
        return '.'.join(sorted(new_reactants_smi)) + '>' + '.'.join(sorted(filter(None, reagents))) + '>' + new_products_smi
    else:
        return None

