# Random walk bridge on trees!
The motivation is as follows: It's clear that a simple random walk on the free product G = Z2*Z2*...*Z2
is the non-negative random walk bridge on the integers. In particular, a bridge
of duration n reaches distance roughly sqrt(n), and "hangs out" there before returning at the end.

I was able to show that for any random walk on G, the bridge either reaches distance sqrt(n) or log(n).
Clearly, the answer should be sqrt(n), but to verify this we run this simulation.

A random walk on a tree returns to the origin with probability exponentially small in the number of
steps. Hence, we can't sample greedily.

We can compute the exact transition probabilities recursively, and sample exactly from the path,
but this would require a stretched exponential amount of space in the number of steps, and
I don't know how to control floating point error. (Assuming rational transition probabilities would
allow for exact calculations using integers, but the denominators would grow like crazy!)

So we implement a random walk bridge by encoding it using a labelled, ordered tree,
where the labels represent generators of G. We define a Markov chain on the labelled, ordered tree
whose stationary measure encodes a uniform word representing the identity in G, ie. a simple random
walk bridge. By applying a Metropolis-Hastings check, we can alter the stationary measure to be that
of any nearest neighbour random walk bridge.

The structure encodes a random walk bridge as follows:
- Each node has a label, which corresponds to a generator of G.
- Perform a depth-first traverse on the tree, and record the labels ak on the nodes in order.
- The process a1*a2*a3***ak thus corresponds to a nearest neighbour walk on G.
- The tree-structure of the Cayley graph of G, along with the fact that each generator of G is
  of order 2, ensures that the process returns to the identity element at the end of the traversal.


# TODO
When I wrote this, I didn't know about Gibbs sampling. Metropolis-Hastings sometimes slows down
a lot when the rejection rate is high. Gibb's sampling would take care of this!
