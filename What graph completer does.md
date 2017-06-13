# The algorithm as it appears in the code

This is the order in which the major algorithm operates.

## 1. Deal with Antipodal stuff.

The antipodal case is strange and great care needs to be taken with it. Adding an edge between different delta edges means we're forced to add another 'relfected' edge.

### 1.1. Make antipodally symmetric if needed.

Make sure every node has an antipode, then add all 'reflected' edges to make all antipodal quadruples symmetric. Pairs of delta edges will end up being completed to a full K4, or will be left as just the two edges for now.

### 1.2. Complete one pode with reduced parameters.

We look at one antipode from each pair and complete this pode using smaller parameters (and other parts of the algorithm). In this case the smaller parameters won't be antipodal.

### 1.3. Make antipodally symmetric if needed.

Extend the edges given in 1.2 to full antipodal quadruples. (This only happens if edges were added in step 1.2.)

## 2. Deal with bipartite stuff.

The bipartite case has a strange way to connect disjoint components.

### 2.1. Complete individual components.

This is just step 3 applied to individual connected components.

### 2.2. Connect components.

Pick a representative from each connected component and add M or M+1 depending on some parity conditions of paths to the representatives.

## 3. Complete forks as perscribed.

The bulk of the algorithm is in completing forks at this stage in a prescribed order.

## 4. Add magic distance (+- 1) to finish.

Once the forks have been dealt with, M or M+1 is added to all remaining non-edges.

### 4.1. Print step 4.

# What happends most of the time

The 3-constrained case (i.e. non-bipartite, non-antipodal) is the most common. Here only steps 3 and 4 are run.
