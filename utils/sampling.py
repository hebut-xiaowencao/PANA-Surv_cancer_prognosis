import numpy as np


class KFold:
    def __init__(self, n_splits=5, shuffle=False, random_state=None):
        if n_splits < 2:
            raise ValueError("n_splits must be at least 2")
        self.n_splits = int(n_splits)
        self.shuffle = shuffle
        self.random_state = random_state

    def get_n_splits(self, X=None):
        return self.n_splits

    def split(self, X):
        n_samples = len(X)
        if self.n_splits > n_samples:
            raise ValueError(f"n_splits={self.n_splits} cannot be greater than n_samples={n_samples}")

        indices = np.arange(n_samples)
        if self.shuffle:
            rng = np.random.default_rng(self.random_state)
            rng.shuffle(indices)

        fold_sizes = np.full(self.n_splits, n_samples // self.n_splits, dtype=int)
        fold_sizes[: n_samples % self.n_splits] += 1

        current = 0
        for fold_size in fold_sizes:
            start, stop = current, current + fold_size
            test_index = indices[start:stop]
            train_index = np.concatenate([indices[:start], indices[stop:]])
            yield train_index, test_index
            current = stop


def resample(values, n_samples, replace=False, random_state=None):
    values = np.asarray(values)
    rng = np.random.default_rng(random_state)
    indices = rng.choice(len(values), size=n_samples, replace=replace)
    return values[indices]
