from abc import abstractmethod

class Splitter:

    def __init__(self, cv: int):
        self.cv = cv
        self.cv_split_id = 0
        self.folds = None

    @abstractmethod
    def split(self):
        pass

    def __iter__(self):
        self.cv_split_id = 0
        return self

    def __next__(self):

        if self.cv_split_id >= self.limit:
            raise StopIteration

        valid_fold = self.folds[self.cv_split_id]
        valid_molecules_id = valid_fold["molecules_id"]
        train_molecules_id = []

        for fold_id, fold_data in self.folds.items():
            if fold_id == self.cv_split_id:
                next

            train_molecules_id = train_molecules_id + fold_data["molecules_id"]

        self.cv_split_id += 1

        return self.cv_split_id, train_molecules_id, valid_molecules_id