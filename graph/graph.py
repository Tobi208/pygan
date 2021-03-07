class Graph:

    def __init__(self):

        self.root = None
        self.nodes = []
        self.edges = []

    def get_node(self, node_id):
        for v in self.nodes:
            if v.node_id == node_id:
                return v