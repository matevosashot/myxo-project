

class SparseTensor:

    def __init__(self, index_map, values):
        self.index_map = index_map
        self.values = values

    def __getitem__(self, index):
        return self.values[self.index_map[index]]
