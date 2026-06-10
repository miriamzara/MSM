from sklearn.decomposition import PCA
class PCAWrapper:
    """
    Small wrapper that gives sklearn PCA a deeptime-like API:
    fit(...).fetch_model().transform(...).
    """
    def __init__(self, dim: int):
        self.model = PCA(n_components=dim)

    def fit(self, traj):
        self.model.fit(traj)
        return self

    def fetch_model(self):
        return self.model

    def transform(self, X):
        return self.model.transform(X)