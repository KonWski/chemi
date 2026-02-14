from .splitter import Splitter
import random

class RandomSplitter(Splitter):
    
    def __init__(self, cv: int, seed=123):
        super().__init__(cv)
        self.seed = seed

    def split(self, molecule_ids, smiles, labels):

        n_obs = len(smiles)
        n_obs_cv = int(n_obs / self.cv)

        random.seed(self.seed)
        internal_molecule_ids = [id for id in range(n_obs)]
        random.shuffle(internal_molecule_ids)

        # folds as fold_id, smiles, labels, original_id
        folds = {fold_id: {"smiles": [], "labels": [], "molecules_id": []} for fold_id in range(self.cv)}

        for fold_id in range(self.cv):

            index_start = fold_id * n_obs_cv
            index_end = (fold_id + 1) * n_obs_cv
            select_molecules_id = internal_molecule_ids[index_start: index_end]

            folds[fold_id]["smiles"] = [smiles[molecule_id] for molecule_id in select_molecules_id]
            folds[fold_id]["labels"] = [labels[molecule_id] for molecule_id in select_molecules_id]
            folds[fold_id]["molecule_ids"] = [molecule_ids[molecule_id] for molecule_id in select_molecules_id]

        self.folds = folds