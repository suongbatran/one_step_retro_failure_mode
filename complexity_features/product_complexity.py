import sys
import os
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from openbabel import openbabel
openbabel.obErrorLog.SetOutputLevel(0)
from complexity_features.boettcher import BottchScorer

from rdkit.Chem import SpacialScore

sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))
from SA_Score import sascorer
from NP_Score import npscorer
fscore = npscorer.readNPModel()

def get_NPScore(mol):
    try:
        score = npscorer.scoreMol(mol,fscore)
        return score
    except Exception:
        return None

def get_SAScore(mol):
    try:
        score = sascorer.calculateScore(mol)
        return score
    except Exception:
        return None

def get_SPScore(mol):
    try:
        score = SpacialScore.SPS(mol, normalize = True)
        return score
    except Exception:
        return None

def get_BoettcherScore(smiles):
    """
    Adapted from https://gitlab.com/mlpds_mit/askcosv2/molecular_complexity.git

    return None if error
    """

    obmol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    obConversion.ReadString(obmol, smiles)

    scorer = BottchScorer(obmol, verbose=False)
    score = scorer.score(obmol)
    
    return score

# smi = "[NH2:1][CH2:2][CH2:3][c:4]1[c:5]2[n:6]([c:7]3[cH:8][cH:9][cH:10][cH:11][c:12]13)[C:13](=[O:14])[CH2:15][CH2:16][CH2:17]2"
# print(get_BoettcherScore(smi))
# print(get_NPScore(Chem.MolFromSmiles(smi)))
# print(get_SAScore(Chem.MolFromSmiles(smi)))
# print(get_SPScore(Chem.MolFromSmiles(smi)))