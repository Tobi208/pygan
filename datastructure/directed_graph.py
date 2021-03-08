class DirectedGraph:
    """
    Represents a directed graph with directed nodes/edges
    """
    def __init__(self):
        self.root = None
        self.nodes = []
        self.edges = []

    def get_node(self, node_id):
        """
        Get node by id

        :param node_id: id of node
        :return: node with id
        """
        for n in self.nodes:
            if n.node_id == node_id:
                return n


class DirectedNode:
    """
    Represents a directed node with ingoing/outgoing edges
    """
    def __init__(self, node_id):
        self.node_id = node_id
        self.in_edges = []
        self.out_edges = []


class DirectedEdge:
    """
    Represents a directed edge with source/target directed nodes
    """
    def __init__(self, source, target):
        self.source = source
        self.target = target