#!/usr/bin/env python
# coding: utf-8

# In[30]:


from itertools import combinations

load("derivative.sage")
load("matroidzeta.sage")


# In[31]:


# Relax a matroid
# This is honestly terrible code; it switches between containers a billion times
def relax_circuit_hyperplane(M):

    # Get circuits and hyperplanes
    circuits = [tuple(c) for c in sorted([sorted(C) for C in M.circuits()])]
    hyperplanes = [tuple(h) for h in sorted([sorted(F) for F in M.hyperplanes()])]
    bases = [sorted (B) for B in M.bases()]

    # Find circuit-hyperplanes
    circuit_hyperplanes = list(set(circuits).intersection(hyperplanes))

    # Return if nothing to relax
    if not circuit_hyperplanes:
        print('No circuit hyperplanes to relax!')
        return None
    
    # Print circuit-hyperplanes
    for ch in enumerate(circuit_hyperplanes):
        print('{0}: '.format(ch[0]) + str(ch[1]))

    # Get circuit-hyperplane to relax
    ch = list(circuit_hyperplanes[int(input())])
    
    #print(bases)

    # Add it as a basis and create new list
    bases.append(ch)
    #new_bases = [''.join(l) for l in bases]
    
    return Matroid(bases=bases)

# cdgh is the circuit-hyperplane to relax to get Vamos from NonVamos
# NV = matroids.named_matroids.NonVamos()
# V = relax_circuit_hyperplane(NV)
# if V.is_isomorphic(matroids.named_matroids.Vamos()):
#    print('NonVamos successfully relaxed to Vamos')

# Get Pappus
# P = matroids.named_matroids.Pappus()

# Plot Pappus to find to relax
# pos = {'a':(0,2), 'b': (1,2), 'c':(2,2), 'd':(0,0), 'e':(1,0), 'f':(2,0), 'g':(.5,1),'h':(1,1), 'i':(1.5,1)}
# P.plot(pos_method=1, pos_dict=pos)

# Relax Pappus and check
# NP = relax_circuit_hyperplane(P)
# if NP.is_isomorphic(matroids.named_matroids.NonPappus()):
#    print('Pappus successfully relaxed to NonPappus')


# In[32]:


def dependent_sets_zeta(M, zeta=False):
    
    # Print relevant information
    print('Dependent sets per cardinality for M({0},{1}):'.format(M.rank(), M.size()), 
          [len(M.dependent_r_sets(r)) for r in range(0, M.rank()+1)])
    
    if zeta:
        print('TZF for M: ', tzf(M))


# In[33]:


g = Matroid(Graph([(0,1,1), (0,1,2), (0,1,3), (1,2,1), (2,3,1)], multiedges=True))
dependent_sets_zeta(g, True)
compare_derivatives(3, g)


# In[42]:


# Gets all matroids of a given rank and ground set size (up to isomorphism)
# by checking all possibe basis combinations 
def get_matroids(r, n):
    
    res = []
    res_2 = []
    
    # from itertools
    sets = combinations(list(range(n)), r)
    
    # Iterate through power set of subsets
    for b in powerset(sets):
        
        # Skip empty
        if not b:
            continue
    
        # Make matroid, see if it is valid
        M = Matroid(bases=b)
        if M.is_valid() and M.size() == 5:
            res.append(M)
            
    # Remove isomorphic copies
    while len(res) > 0:
        res_2.append(res[0])
        res = [m for m in res[1:] if not m.is_isomorphic(res[0])]

    return res_2

# Get the set of dependent matroid we would want
M = get_matroids(3,5)
M2 = [m for m in M if len(m.dependent_r_sets(2)) > 0]
  
for i, m in enumerate(M2):
    print('Matroid index: {0}'.format(i))
    compare_derivatives(3, m)

