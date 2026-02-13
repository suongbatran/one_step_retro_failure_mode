from rxnmapper import BatchedMapper
rxn_mapper = BatchedMapper(batch_size=256)
# rxns = ['CC[O-]~[Na+].BrCC>>CCOCC', 'invalid>>reaction', 'C1COCCO1.CC(C)(C)[O:3][C:2](=[O:1])[CH2:4][O:5][NH:6][C:7](=[O:8])[NH:9][CH2:10][c:11]1[cH:12][cH:13][cH:14][c:15]2[cH:16][cH:17][cH:18][cH:19][c:20]12.Cl>>O=C(O)CONC(=O)NCc1cccc2ccccc12']
rxns =["C1CCOC1.NCc1cccnc1.O=C(O)C1Oc2ccc(Cl)cc2Cc2cc(Cl)ccc2O1>>O=C(NCc1cccnc1)C1Oc2ccc(Cl)cc2Cc2cc(Cl)ccc2O1"]
# The following calls work with input of arbitrary size. Also, they do not raise 
# any exceptions but will return ">>" or an empty dictionary for the second reaction.
results = list(rxn_mapper.map_reactions(rxns))  # results as strings directly
print(results)

results = list(rxn_mapper.map_reactions_with_info(rxns))  # results as dictionaries (as above)
print(results)