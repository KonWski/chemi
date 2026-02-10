from .splitter import Splitter
import random

class ScaffoldSplitter(Splitter):
    
    def __init__(self, cv: int, seed):
        super().__init__(cv)
        self.seed = seed

    def split(self, smiles, labels):

        random.seed(self.seed)
        molecules_id = [molecule_id for molecule_id in range(len(smiles))]
        shuffled_molecules_id = random.shuffle(molecules_id)

        n_obs = len(molecules_id)
        n_obs_cv = int(n_obs / self.cv)

        # folds as fold_id, smiles, labels, original_id
        folds = {fold_id: {"smiles": [], "labels": [], "molecules_id": []} for fold_id in range(self.cv)}

        for fold_id in range(self.cv):

            index_start = fold_id * n_obs_cv
            index_end = (fold_id + 1) * n_obs_cv
            select_molecules_id = shuffled_molecules_id[index_start: index_end]

            folds[fold_id]["smiles"] = [smiles[molecule_id] for molecule_id in select_molecules_id]
            folds[fold_id]["labels"] = [labels[molecule_id] for molecule_id in select_molecules_id]
            folds[fold_id]["molecules_id"] = select_molecules_id

        self.folds = folds