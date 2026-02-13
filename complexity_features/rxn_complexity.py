from rdkit import Chem


def count_num_react_atom(rxn_smiles):
    """
    Count number of reacting atoms in reactants
    """
    
    try:
        parts = rxn_smiles.split('>')
        if len(parts) != 3:
            return -1 
    except:
        return -1

    reactants = Chem.MolFromSmiles(parts[0])
    products = Chem.MolFromSmiles(parts[2])
    if reactants is None or products is None:
        return -1

    # Map atom map numbers to atoms (excluding unmapped atoms)
    reactant_atoms = {a.GetAtomMapNum(): a for a in reactants.GetAtoms() if a.GetAtomMapNum() != 0}
    product_atoms = {a.GetAtomMapNum(): a for a in products.GetAtoms() if a.GetAtomMapNum() != 0}

    atom_changed = 0

    for mapnum in reactant_atoms:
        r_atom = reactant_atoms[mapnum]
        p_atom = product_atoms[mapnum]

        r_types = sorted(nbr.GetSymbol() for nbr in r_atom.GetNeighbors())
        p_types = sorted(nbr.GetSymbol() for nbr in p_atom.GetNeighbors())

        r_neighbors = sorted(n.GetAtomMapNum() for n in r_atom.GetNeighbors())
        p_neighbors = sorted(n.GetAtomMapNum() for n in p_atom.GetNeighbors())

        # 1. Compare formal charge
        if r_atom.GetFormalCharge() != p_atom.GetFormalCharge():
            atom_changed += 1
            # also count unmapped neighbors
            atom_changed += r_neighbors.count(0)
            continue

        # 2. Compare neighbors' atom type

        if r_types != p_types:
            atom_changed += 1
            atom_changed += r_neighbors.count(0)
            # print(atom_changed)
            continue

        # 3. Compare neighbors' atom map numbers
        if r_neighbors != p_neighbors:
            atom_changed += 1
            atom_changed += r_neighbors.count(0) 
            continue

        # 4. Compare bond types to neighbors
        r_bonds = sorted(
            (n.GetAtomMapNum(), reactants.GetBondBetweenAtoms(r_atom.GetIdx(), n.GetIdx()).GetBondType())
            for n in r_atom.GetNeighbors()
        )
        p_bonds = sorted(
            (n.GetAtomMapNum(), products.GetBondBetweenAtoms(p_atom.GetIdx(), n.GetIdx()).GetBondType())
            for n in p_atom.GetNeighbors()
        )
        if r_bonds != p_bonds:
            atom_changed += 1

    return atom_changed

def count_num_react_atom_retro(rxn_smiles):
    """
    Count number of reacting atoms in products
    """
    
    try:
        parts = rxn_smiles.split('>')
        if len(parts) != 3:
            return -1  # Invalid reaction SMILES format
    except:
        return -1

    reactants = Chem.MolFromSmiles(parts[0])
    products = Chem.MolFromSmiles(parts[2])
    if reactants is None or products is None:
        return -1

    # Map atom map numbers to atoms (excluding unmapped atoms)
    reactant_atoms = {a.GetAtomMapNum(): a for a in reactants.GetAtoms() if a.GetAtomMapNum() != 0}
    product_atoms = {a.GetAtomMapNum(): a for a in products.GetAtoms() if a.GetAtomMapNum() != 0}


    atom_changed_list = []

    for mapnum in product_atoms:
        if mapnum not in reactant_atoms:
            atom_changed_list.append(mapnum)
            # mapped numbers in products but not in reactants are considered as reaction centers 
            # (there should be this cases if use "matching" mapping style)
            continue 
        
        r_atom = reactant_atoms[mapnum]
        p_atom = product_atoms[mapnum]


        r_types = sorted(nbr.GetSymbol() for nbr in r_atom.GetNeighbors())
        p_types = sorted(nbr.GetSymbol() for nbr in p_atom.GetNeighbors())

        r_neighbors = sorted(n.GetAtomMapNum() for n in r_atom.GetNeighbors())
        p_neighbors = sorted(n.GetAtomMapNum() for n in p_atom.GetNeighbors())

        # 1. Compare formal charge
        if p_atom.GetFormalCharge() != r_atom.GetFormalCharge():
            atom_changed_list.append(mapnum)
            continue

        # 2. Compare neighbors' atom type
        if p_types != r_types:
            atom_changed_list.append(mapnum)
            continue

        # 3. Compare neighbors' atom map numbers
        if p_neighbors != r_neighbors:
            atom_changed_list.append(mapnum)
            continue

        # 4. Compare bond types to neighbors
        r_bonds = sorted(
            (n.GetAtomMapNum(), reactants.GetBondBetweenAtoms(r_atom.GetIdx(), n.GetIdx()).GetBondType())
            for n in r_atom.GetNeighbors()
        )
        p_bonds = sorted(
            (n.GetAtomMapNum(), products.GetBondBetweenAtoms(p_atom.GetIdx(), n.GetIdx()).GetBondType())
            for n in p_atom.GetNeighbors()
        )
        if p_bonds != r_bonds:
            atom_changed_list.append(mapnum)
    return len(atom_changed_list) 

def count_num_change_ring(rxn_smiles):
    """
    Count number of changing rings in reactants and products
    """
    try:
        parts = rxn_smiles.split('>')
        if len(parts) != 3:
            return -1  # Invalid reaction SMILES format
    except:
        return -1

    reactants = Chem.MolFromSmiles(parts[0])
    products = Chem.MolFromSmiles(parts[2])
    if reactants is None or products is None:
        return -1
    
    def get_mapped_rings(mols): 
        # because deprotection rxn does not include the protecting group in products, it would be better to use rings containg all mapped atoms
        rings = set()
        ring_info = mols.GetRingInfo()
        mapped_atoms = {atom.GetIdx() for atom in mols.GetAtoms() if atom.GetAtomMapNum() != 0}
        for ring in ring_info.AtomRings():
            for atom_idx in ring:
                if atom_idx not in mapped_atoms:
                    break
            else:
                mapped_ring = frozenset(sorted(atom.GetAtomMapNum() for atom in mols.GetAtoms() if atom.GetIdx() in ring))
                rings.add(mapped_ring)
        return rings

    reactant_rings = get_mapped_rings(reactants)
    product_rings = get_mapped_rings(products)
    
    # Count rings that are different after reaction
    changed_rings = reactant_rings.symmetric_difference(product_rings)
    # print(changed_rings)
    return len(changed_rings)

# rxn_smiles = "Br[CH2:19][CH2:20][CH2:21][CH2:22][O:23][c:24]1[cH:25][c:26]2[c:31]([cH:32][cH:33]1)[NH:30][C:29](=[O:34])[CH2:28][CH2:27]2.[SH:1][c:2]1[n:3][cH:4][cH:5][cH:6][cH:7]1>CS(C)=O.O.O=C([O-])[O-].[K+].[K+]>[S:1]([c:2]1[n:3][cH:4][cH:5][cH:6][cH:7]1)[CH2:19][CH2:20][CH2:21][CH2:22][O:23][c:24]1[cH:25][c:26]2[c:31]([cH:32][cH:33]1)[NH:30][C:29](=[O:34])[CH2:28][CH2:27]2"
# rxn_smiles = 'O=[C:1]([NH:2][c:4]1[c:8]([CH:12]=O)[cH:13][cH:16][c:14]([C:17]([CH3:18])([CH3:19])[CH3:20])[cH:9]1)[C:3]1([CH3:5])[CH2:6][CH2:10][CH2:15][CH2:11][CH2:7]1.[CH3:21][C:22]1([CH3:23])[c:24]2[c:26]([cH:31][c:35]([NH2:37])[c:32]([NH2:36])[cH:27]2)[C:30]([CH3:33])([CH3:34])[C:25]1([CH3:28])[CH3:29]>>[c:1]1([C:3]2([CH3:5])[CH2:6][CH2:10][CH2:15][CH2:11][CH2:7]2)[n:2][c:4]2[c:8]([c:12]3[n:36]1[c:32]1[cH:27][c:24]4[c:26]([cH:31][c:35]1[n:37]3)[C:30]([CH3:33])([CH3:34])[C:25]([CH3:28])([CH3:29])[C:22]4([CH3:21])[CH3:23])[cH:13][cH:16][c:14]([C:17]([CH3:18])([CH3:19])[CH3:20])[cH:9]2'
# print(count_num_react_atom_retro(rxn_smiles))
# print(count_num_change_ring(rxn_smiles))