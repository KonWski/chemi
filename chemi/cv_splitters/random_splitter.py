from .splitter import Splitter
import random
import numpy as np

class RandomSplitter(Splitter):
    
    def __init__(self, cv: int, seed=123, stratified=True):
        super().__init__(cv)
        self.seed = seed
        self.stratified = stratified


    def split(self, molecule_ids, smiles, labels):

        if self.stratified:
            self.stratified_split(molecule_ids, smiles, labels)
        else:
            self.unstratified_split(molecule_ids, smiles, labels)


    def unstratified_split(self, molecule_ids, smiles, labels):

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
    

    def stratified_split(self, molecule_ids, smiles, labels):

        n_obs = len(smiles)
        internal_molecule_ids = [id for id in range(n_obs)]
        minority_indices = [molecule_id for molecule_id in internal_molecule_ids if labels[molecule_id] == 1]
        majority_indices = [molecule_id for molecule_id in internal_molecule_ids if labels[molecule_id] == 0]

        random.shuffle(minority_indices)
        random.shuffle(majority_indices)

        # folds as fold_id, smiles, labels, original_id
        folds = {fold_id: {"smiles": [], "labels": [], "molecules_id": []} for fold_id in range(self.cv)}
        minority_folds = np.array_split(minority_indices, self.cv)
        majority_folds = np.array_split(majority_indices, self.cv)

        for fold_id in range(self.cv):

            minority_ids = minority_folds[fold_id]
            majority_ids = majority_folds[fold_id]
            fold_molecule_ids = list(minority_ids) + list(majority_ids)

            folds[fold_id]["smiles"] = [smiles[molecule_id] for molecule_id in fold_molecule_ids]
            folds[fold_id]["labels"] = [labels[molecule_id] for molecule_id in fold_molecule_ids]
            folds[fold_id]["molecule_ids"] = [molecule_ids[molecule_id] for molecule_id in fold_molecule_ids]

        self.folds = folds