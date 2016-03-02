import random
import tree
import pdb
from IPython import embed

class Bridge(tree.Tree):
    '''
    Implement a random walk bridge by storing encoding it using a labelled, ordered tree,
    where the labels represent generators from the free product G = Z2*Z2*...*Z2.

    See the readme for explanation of this encodes a random walk bridge.
    '''

    class Node(tree.Tree.Node):

        def __init__(self, data={}):
            super(Bridge.Node, self).__init__(data)
            self.generator = None
            self.candidate_generator = None
            self.label = None
            self.indices_in_parent_list = []

    def __init__(self,
                 num_nodes,
                 threshold=0.0,
                 excursion=False,
                 metropolis=False,
                 generators=['a', 'b', 'c', 'd'],
                 gen_weights=[.25, .25, .25, .25],
                 verbose=False
                 ):
        '''
        excursion:  Determines restriction on the degree of the root
        threshold:  (small) probability that the root's label changes
        metropolis: Use the metropolis algorithm to engineer the stationary measure
        verbose:    Determines whether to print information to stdout.
        '''
        super(Bridge, self).__init__()
        self.generators = generators
        self.gen_weights = gen_weights
        self.weights = dict(p for p in zip(self.generators, gen_weights))

        self.threshold = threshold

        self.root = Bridge.Node()
        self.root.label = 0
        self.root.generator = random.choice(self.generators)
        self.num_generators = len(generators)

        # Create a fake parent of the root. See readme for explanation.
        self.parent_root = Bridge.Node()
        self.root.parent = self.parent_root
        self.parent_root.label = -1
        self.parent_root.generator = None

        # Determine the behaviour of the chain.
        self.metropolis = metropolis
        self.excursion = excursion
        self.verbose = verbose

        if excursion:
            self.root.generator = random.choice(self.generators)
        else:
            self.root.generator = None

        self.nodes = [self.root]
        self.node_indices = range(len(self.nodes))
        self.leaves = {self.root: False}
        self.num_leaves = 0

        # Hard code a limit for now. Optimize later.
        self.num_nodes = num_nodes
        self.normalizing_constant = 2*self.num_nodes - 1


        # We initialize the tree to be one root with many leaves.

        # We store a list of potential parents for a new leaf. This
        # allows us to sample a new position for a leaf in O(1) time,
        # using O(num_nodes) space.
        # To update the list of potential parents in O(1) time, we store
        # the indices of the node in the parent list
        self.parent_list = [self.root for i in range(num_nodes + 1)]
        self.root.indices_in_parent_list = range(num_nodes + 1)

        for i in range(1, self.num_nodes):
            self.add_child(self.root, {'label':i,
                                       'indices':[num_nodes + i + 1]})

        for child in self.root.children:
            self.parent_list.append(child)

        # Store the letter frequencies in the encoded G-word.
        self.word_profile = {}
        # self.update_word_profile()

        self.num_rejections = 0

        # Check our assumptions
        max_num_nodes = 10000
        assert(num_nodes < max_num_nodes)
        assert(self.num_nodes == len(self.nodes))
        assert(sum(self.weights[g] for g in self.weights) == 1)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def add_child(self, parent, data={}):
        '''Adds a new child below the parent node.'''

        child = Bridge.Node()
        child.parent = parent
        parent.children.append(child)
        child.generator = self.get_new_generator(self.root.generator)
        try:
            child.label = data['label']
            child.indices_in_parent_list = data['indices']
            self.nodes.append(child)
            self.leaves[child] = True
            self.num_leaves += 1
        except KeyError:
            print("Please pass in the correct data.")
            exit()

    def get_random_leaf(self):
        '''Chooses a uniform leaf.

        Since the tree looks like a critical Galton Watson tree, the number
        of leaves at stationary will be roughly half the number of nodes, so
        this takes Geom(2) time, or O(1) time in expectation.
        '''
        leaf = random.choice(self.nodes)
        while not self.leaves[leaf]:
            leaf = random.choice(self.nodes)
        return leaf

    def move_random_leaf(self):
        '''
        Moves a randomly chosen leaf to a uniform position. Chooses
        a uniform labelling for the new leaf, satisfying the labelling
        constraints.

        If the metropolis flag is set to True, the move is rejected with
        a certain probability, allowing the specification of the stationary
        measure.
        '''

        # Choose a uniform node to move.
        old_leaf = self.get_random_leaf()
        old_parent = old_leaf.parent

        # Choose a random parent for the new node. The new node is in
        # a uniform position, so the parent should be chosen with
        # probability proportional to len(parent.children) + 1
        new_parent = random.choice(self.parent_list)
        while new_parent == old_leaf:
            new_parent = random.choice(self.parent_list)

        if self.verbose:
            print ('Moving node %i to underneath node %i'
                   % (old_leaf.data['label'], new_parent.data['label']))

        # Find a new relabelling of the tree below the old node. The
        # new relabelling must be consistent with the new parent.
        top = True
        self.candidate_relabel(old_leaf, new_parent.generator, top)

        if old_leaf.candidate_generator == new_parent.generator:
            raise ValueError('candidate_generator is broken')

        # Relabel the tree under the new node.
        if self.metropolis:
            '''
            Calculate the transition threshold.
            Note that if we are sampling an excursion, and not a
            bridge, that the markov chain is in fact symmetric! So,
            we don't need to worry about anything other than the
            change in the stationary measure.
            '''

            if self.verbose:
                print '''The candidate change:
                    old leaf: \t\t%i
                    old parent: \t%i
                    new parent: \t%i
                    new parent's generator: %s
                    old generator:          %s
                    candidate generator:    %s
                ''' % (old_leaf.data['label'],
                       old_leaf.parent.data['label'],
                       new_parent.data['label'],
                       new_parent.generator,
                       old_leaf.generator,
                       old_leaf.candidate_generator)

            acceptance_prob = 1.0
            p = self.weights[old_leaf.candidate_generator] ** 2
            q = self.weights[old_leaf.generator] ** 2
            acceptance_prob *= p/q

            # Adjust for the edge case when exactly one of the
            # parents is the root.

            #TODO check this edge case.
            if not self.excursion:
                d = self.num_generators
                if old_parent == self.root and new_parent != self.root:
                    acceptance_prob *= 1.0*d/(d-1)
                elif new_parent == self.root and old_parent != self.root:
                    acceptance_prob *= 1.0*(d-1)/d

        else:
            acceptance_prob = 1.0

        acceptance_prob = min(acceptance_prob, 1.0)

        x = random.random()
        if x < acceptance_prob:
            self.relabel(old_leaf)
            self.accept_move(old_leaf, new_parent)

        else:
            self.num_rejections += 1
            if self.verbose:
                print ('Rejected. x = %f, acceptance_prob = %f'
                       % (x, acceptance_prob))

    def accept_move(self, old_leaf, new_parent):

        # Update parent_list.
        old_parent_index = old_leaf.parent.indices_in_parent_list[-1]
        del old_leaf.parent.indices_in_parent_list[-1]

        self.parent_list[old_parent_index] = new_parent
        new_parent.indices_in_parent_list.append(old_parent_index)
        
        # Delete old node.
        old_leaf.parent.children.remove(old_leaf)
        if len(old_leaf.parent.children) == 0:
            self.leaves[old_leaf.parent] = True

        # Insert new node.
        old_leaf.parent = new_parent
        self.leaves[old_leaf.parent] = False
        index = random.randrange(len(new_parent.children)+1)
        new_parent.children.insert(index, old_leaf)

        if self.verbose:
            print 'The new tree structure: '
            tree.verbose_traverse()

            print 'Node %i has been moved.' % old_leaf.label
            print '-'*20

    def print_labels(self, list_of_nodes):
        print(map(lambda node: node.label, list_of_nodes) )


    def run(self, num_iterations):
        for i in range(num_iterations):
            # if i % (num_steps/100) == 0 and i > 0:
                # print '%i percent complete.' % (i*100/num_steps)
            if self.verbose:
                print 'Next step:'
            x = random.random()
            if x < self.threshold:
                # Relabel the root randoml
                if self.verbose:
                    print 'Changing the root label.'
                top = True
                self.candidate_relabel(self.root, None, top)

                threshold = self.get_label_threshold()

                x = random.random()
                if x < threshold:
                    self.relabel(self.root)

                if self.verbose:
                    print 'The root label has changed! The tree structure is:'
                    self.verbose_traverse()
            else:
                # Move a random leaf
                self.move_random_leaf()
            if self.verbose:
                if i % 100 == 1:
                    print 'The tree structer after %i steps is: ' % i
                    self.verbose_traverse()

                print 'Checking the labels.'
                if not self.check_labels():
                    print 'Check the labels:'
                    self.verbose_traverse()
                print 'Labels are all good.'

                print 'Step complete.'
                print '-'*20

    def candidate_relabel(self, node, invalid_generator=None, top=False):
        '''Updates the generator at a node with a randomly chosen
        generator. relabel is called recursively on all children,
        resulting in a (candidate) uniform random labeling of the tree
        below node.
        '''

        # Complexity: unclear. Should be O(sqrt n) on the average

        # Different behaviour when node is the root:
        if node == self.root:
            node.candidate_generator = random.choice(self.generators)
            for child in node.children:
                self.candidate_relabel(child, node.candidate_generator)
            return

        node.candidate_generator = self.get_new_generator(invalid_generator)

        if self.verbose:
            print ('''The candidate relabellings for node %i and its parent,
                   node %i, is:  %s and %s'''
                   % (node.label,
                      node.parent.label,
                      node.candidate_generator,
                      node.parent.candidate_generator
                      )
                   )

        for child in node.children:
            self.candidate_relabel(child, node.candidate_generator)

    def relabel(self, node):
        '''
        Sets node.generator to be node.candidate_generator. Relabel is
        called recursively on all of the children of node, relabelling
        all of the nodes below to the new labelling.
        '''

        if self.verbose:
            print 'Relabeling node \t%i.' % node.label
            print 'old generator: \t%s' % node.generator
            print 'new generator: \t%s' % node.candidate_generator
            if node != self.root:
                print 'parent\'s gen: \t%s' % node.parent.generator
        node.generator = node.candidate_generator

        for child in node.children:
            self.relabel(child)

    def get_label_threshold(self):
        threshold = 1.0
        for node in self.nodes:
            p = self.weights[node.candidate_generator] ** 2
            q = self.weights[node.generator] ** 2
            threshold *= p/q
        return threshold

    def verbose_traverse(self):
        traversal = self.traverse(self.root)

        for node in traversal:
            if node != self.root:
                print node.label, node.generator
            else:
                print node.label, 'root'

    def tree_to_word(self):
        '''Converts the labelled tree into the word in the free group
        represented by the tree.'''
        word = ''

        traversal = self.traverse(self.root)

        old_node = traversal.next()
        word = self.root.generator
        new_node = traversal.next()
        while True:
            if new_node.generator and new_node != old_node.parent:
                word = word + new_node.generator
            if old_node.parent == new_node and old_node.generator:
                word = word + old_node.generator
            # if not new_node.children:
                word = word + new_node.generator
            try:
                old_node, new_node = new_node, traversal.next()
            except StopIteration:
                break
        word = word + self.root.generator
        return word

    def check_normalizing_constant(self):
        total = 0
        for node in self.nodes:
            total += len(node.children) + 1

        # print total, self.normalizing_constant

    def check_labels(self):
        for node in self.nodes[1:]:
            if node.generator == node.parent.generator:
                print ('''Node %i has the same label as node %i\'s\
                       parent: ''' % (node.label, node.label)
                       )
                print node.generator, node.parent.generator
                raise ValueError('candidate_relabel is broken.')
                return False
        return True

    def get_new_parent(self, node):
        '''
        Finds a random parent of node. Parent is chosen with probability
        proportional to len(parent.children) + 1, so that inserting a
        new node in a uniform position below parent is equivalent to
        selecting a uniform position for a new leaf in the tree.
        '''

        new_parent = random.choice(self.parent_list)

        while not self.valid_new_parent(new_parent, node):
            new_parent = random.choice(self.parent_list)

        return new_parent

    def valid_new_parent(self, parent, child):
        '''Checks if parent is in the subtree below child.'''

        traversal = self.traverse(child)
        for node in traversal:
            if node == parent:
                return False

        return True

    def update_word_profile(self):
        self.word_profile = dict((g, 0) for g in self.generators)
        total_letters = 0
        for node in self.nodes[1:]:
            self.word_profile[node.generator] += 1
            total_letters += 1
        for g in self.generators:
            self.word_profile[g] *= 1.0/total_letters

    def get_new_generator(self, invalid_generator, max_tries=100000000000):
        new_generator = random.choice(self.generators)

        # This is inexplicably an infinite loop without the max_tries
        # counter. However, the running time does not seem to depend on
        # the actual value of max_tries.
        while new_generator == invalid_generator and max_tries:
            new_generator = random.choice(self.generators)
            max_tries = max_tries - 1

        return new_generator


def sample_from_chain(generators,
                      gen_weights,
                      num_nodes,
                      num_steps,
                      threshold,
                      excursion,
                      metropolis,
                      verbose
                      ):
    '''Samples from the chain. Prints the output.'''

    import time
    time1 = time.time()

    print 'Generators: ', generators
    print 'Weights: ', gen_weights

    import math
    # num_steps = int(1000*math.log(num_nodes)*num_nodes)
    # num_steps = min(num_steps, 2000)
    # num_steps = int(raw_input('Number of steps: '))
    # num_steps = 100

    print 'Threshold (for changing the root label): %f' % threshold

    tree = Bridge(num_nodes,
                  threshold,
                  excursion,
                  metropolis,
                  generators,
                  gen_weights,
                  verbose)

    # print 'After %i steps:' % num_steps
    # tree.verbose_traverse()
    #
    #
    # num_steps = 200
    # for i in range(num_steps):
    #     tree.run(1)
    #     word = tree.tree_to_word()
    #     if word[0] == word[1]:
    #         print 'problem: '
    #         tree.verbose_traverse()
    #         print word
    #
    print 'Number of nodes: %i' % num_nodes
    print 'Total number of steps:    %i' % num_steps

    tree.run(num_steps)
    print 'Number of rejected steps: %i' % tree.num_rejections
    print 'Number of accepted steps: %i' % (num_steps - tree.num_rejections)

    height = tree.update_height(True)

    normalized_height = 1.0*height/math.sqrt(num_nodes)
    print ('The normalized height of the tree is: %f * sqrt(num_steps)'
           % normalized_height)

    print 'The root has %i children. ' % len(tree.root.children)

    # The incorrect profile
    wa, wb, wc = (w*(1-w) for w in gen_weights)
    zw = wa + wb + wc
    wa, wb, wc = wa/zw, wb/zw, wc/zw

    # Probabilities for the relevant non-backgracking random walk
    pa, pb, pc = (g**2 for g in gen_weights)
    zp = pa + pb + pc
    pa, pb, pc = pa/zp, pb/zp, pc/zp

    # The profile for the non-backtracking random walk
    qa, qb, qc = pa*(1-pa), pb*(1-pb), pc*(1-pc)
    zq = qa + qb + qc
    qa, qb, qc = (qa/zq, qb/zq, qc/zq)

    tree.update_word_profile()
    print 'The word profile is:'
    print 'Expected\t\t\tActual\t\t\tWrong'
    print qa, '\t\t', tree.word_profile[generators[0]], '\t', wa
    print qb, '\t\t', tree.word_profile[generators[1]], '\t', wb
    print qc, '\t\t', tree.word_profile[generators[2]], '\t', wc

    time2 = time.time()
    print 'Elapsed time: %f seconds' % (time2 - time1)

    # print 'The tree structure is:'
    # tree.verbose_traverse()

    # word = tree.tree_to_word()
    # print '\n\n\n'
    # print 'And finally, the random word is: %s' % word

if __name__ == '__main__':

    generators = ['A', 'B', 'C']
    # gen_weights = [.2, .2, .3, .3]
    # gen_weights = [.25, .25, .25, .25]
    # gen_weights = [1.0/3, 1.0/3, 1.0/3]
    gen_weights = [.2, .2, .6]

    num_nodes = 5000
    num_steps = 500

    threshold = 0.000001
    # threshold = 1000.0/num_steps
    excursion = True
    metropolis = True
    verbose = False

    params = ('(%r, %r, %i, %i, %f, %r, %r, %r)' %
             (generators,
              gen_weights,
              num_nodes,
              num_steps,
              threshold,
              excursion,
              metropolis,
              verbose
              )
              )

    def check(b):
        t = b.traverse(b.root)
        count = {}
        for n in t:
            count[n] = len(n.indices_in_parent_list)

        print sum(count[n] for n in count)

    b = Bridge(num_nodes,
               threshold,
               excursion,
               metropolis,
               generators,
               gen_weights,
               verbose
               )

    embed()

    # import profile
    # profile.run('sample_from_chain%s' % params)
