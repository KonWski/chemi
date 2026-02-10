from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from collections import Counter
import heapq
import numpy as np
from .splitter import Splitter

class ScaffoldSplitter(Splitter):
    
    def __init__(self, cv: int):
        super().__init__(cv)

    def split(self, molecule_ids, smiles, labels):
        # TODO molecules_id odpowiadaja wewnetrznemu oznaczeniu a nie idkom z pierwotnego datasetu

        scaffolds = []
        for smile in smiles:
            mol = Chem.MolFromSmiles(smile)
            mol = Chem.AddHs(mol)
            scaffolds.append(Chem.MolToSmiles(GetScaffoldForMol(mol)))

        # smile for each scaffold id created from unique scaffolds
        indexed_scaffolds = {smile: id for id, smile in enumerate([np.unique(scaffolds)])}
        
        # scaffold id for each molecule
        scaffold_ids = [indexed_scaffolds[smile] for smile in scaffolds]

        # group scaffold into balanced groups
        # unique_scaffold_id -> fold_id
        grouped_scaffold_ids = self.reverse_dict(self.group_scaffolds_balanced(scaffold_ids, self.cv))

        # folds as fold_id, smiles, labels, original_id
        folds = {fold_id: {"smiles": [], "labels": [], "molecules_id": []} for fold_id in range(self.cv)}

        for scaffold_iter, scaffold_id in enumerate(scaffold_ids):
            fold_id = grouped_scaffold_ids[scaffold_id]
           
            folds[fold_id]["smiles"].append(smiles[scaffold_iter])
            folds[fold_id]["labels"].append(labels[scaffold_iter])
            folds[fold_id]["molecule_ids"].append(molecule_ids[scaffold_iter])

        self.folds = folds


    def group_scaffolds_balanced(self, scaffold_ids):
        
        # number of compounds with the same scaffold
        counts = Counter(scaffold_ids)

        # folds as (current_size, bucket_index, [scaffold_ids])
        folds = [(0, fold_index, []) for fold_index in range(self.cv)]
        
        # convert the folds into a heap
        heapq.heapify(folds)

        # sort ids by frequency (largest first)
        for id, freq in sorted(counts.items(), key=lambda x: -x[1]):
            
            # smallest element in the heap
            size, fold_index, scaffold_ids = heapq.heappop(folds)

            # add the smallest scaffolds group to the largest one
            scaffold_ids.append(id)

            # update the heap with extended number of observations
            heapq.heappush(folds, (size + freq, fold_index, scaffold_ids))

        # fold_id: scaffolds id belonging to it
        return {fold_id: set(scaffold_ids) for _, fold_id, scaffold_ids in folds}


    def reverse_dict(self, d):
        "Unpacks lists stored as values and uses them as new keys"

        reversed_dict = {}

        for original_key, list_values in d.items():
            for list_element in list_values:
                reversed_dict[list_element] = original_key

        return reversed_dict