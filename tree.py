
class Tree(object):
    '''
    A rooted tree data structure with some tree manipulation methods.

    Data may be stored at each node. This allows for eg. a coloured tree, or
    labelled vertices.
    '''

    class Node(object):
        """A node in a tree. Contains pointers to the parent, as well as
           an arbitrary number of children."""
        def __init__(self, data=None):
            self.data = data
            self.parent = None
            self.children = []
            self.height = 0

    def __init__(self):
        '''Initialize the tree to contain a single node (the root)'''
        self.root = Tree.Node()

        # maintain a list of all the nodes to speed up random node sampling
        self.nodes = [self.root]

    def add_child(self, parent, child_data):
        '''Add a new child below the parent node.'''

        child = Tree.Node(child_data)
        child.parent = parent
        parent.children.append(child)
        self.nodes.append(child)

    def delete_branch(self, node):
        '''Delete the branch rooted at node.'''
        self.node.parent.children.delete(node)

        # Recursively delete offspring.
        traversal = self.traverse(node, False)

        for node in traversal:
            self.nodes.remove(node)

    def move_child(self, child, new_parent, at_random=False):
        '''Changes child's parent from parent to new_parent.'''

        child.parent.children.remove(child)
        child.parent = new_parent

        if not at_random:
            child.parent.children.append(child)

        if at_random:
            pass

    def get_children_data(self, node):
        for child in node.children:
            yield child.data

    def update_height(self, return_max_height=False):
        '''Traverses the tree, and updates the heights.

        This is not an efficient way to maintain the height of the tree,
        as it performs a full depth-first traversal.'''

        traversal = self.traverse(self.root)

        max_height = 0
        for node in traversal:
            if not node.height and node != self.root:
                node.height = node.parent.height + 1

            if node.height > max_height:
                max_height = node.height

        if return_max_height:
            return max_height

    def traverse(self, current_node, repeats=True):
        '''Perform a depth first traversal, yielding a node each time
           it's visited, starting from current_node.'''
        yield current_node
        for child in current_node.children:
            traverse = self.traverse(child)
            for node in traverse:
                yield node
            yield current_node

    def verbose_traverse(self, current_node):
        '''Prints out the tree structure.

        If the nodes have unique labels, this allows the reconstruction
        of the tree.'''

        traversal = self.traverse(current_node)
        for node in traversal:
            print node.data
