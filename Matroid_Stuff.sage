#!/usr/bin/env python
# coding: utf-8

# In[1]:


from itertools import combinations

load('derivative.sage')
load('matroidzeta.sage')
load('fan.sage')


# In[2]:


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


# In[3]:


def dependent_sets_zeta(M, zeta=False):
    
    print('Matroid Rank: {0}, Matroid Size: {1}.'.format(M.rank(), M.size()))
    print('Entry (i, j) is the number of cardinality j, rank i dependent sets (one indexed):')
    
    # Get array of 0s
    res = [ [0] * (M.rank()+1) for _ in range(M.rank()) ]
    
    # Get dependent sets
    dep_sets = [M.dependent_r_sets(r) for r in range(0, M.rank()+1)]
    
    # Split dependent sets by rank
    for r, dep_set in enumerate(dep_sets):
        for dep in dep_set:
            res[M.rank(dep)][r] += 1
    
    # Print array
    for row in res[1:]:
        print(row[1:])
    
    if zeta:
        print('TZF for M: ', tzf(M))


# In[4]:


# Displays a matroid
def display_matroid(M, d):
    if M.rank() <= 3:
        M.plot().show()
    compare_derivatives(d, M)
    dependent_sets_zeta(M)
    print('')


# In[5]:


# Gets all matroids of a given rank and ground set size (up to isomorphism)
# by checking all possibe basis combinations 
def get_matroids(r, n):
    
    res = []
    res_2 = []
    
    # from itertools
    sets = combinations(list(range(n)), r)
    
    # Iterate through power set of subsets
    for b in powerset(sets):
        
        # Don't use empty
        if not b:
            continue

        # Make matroid, see if it is valid
        M = Matroid(bases=b)
        if M.size() == n and M.is_valid():
            res.append(M)
            
    # Remove isomorphic copies
    while len(res) > 0:
        res_2.append(res[0])
        res = [m for m in res[1:] if not m.is_isomorphic(res[0])]
    return res_2

# Get the set of dependent matroid we would want
M = get_matroids(4,6)

M2 = [m for m in M if len(m.dependent_r_sets(2)) > 0]
  
for i, m in enumerate(M2):
    print('Matroid index: {0}'.format(i))
    display_matroid(m, 4)


# In[6]:


# Finds D(42) and D(32) for fourth derivative
# Create a graphic matroid
G = Graph([(0,1,1), (0,1,2),(1,2,1),(2,3,1), (2,4,1),(3,4,1)], multiedges=True)
G.plot().show()
M = Matroid(G)
display_matroid(M, 4)

G2 = Graph([(0,1,1), (3,4,2),(1,2,1),(2,3,1), (2,4,1),(3,4,1)], multiedges=True)
G2.plot().show()
M2 = Matroid(G2)
display_matroid(M2, 4)

G3 = Graph([(0,1,1), (0,1,2),(1,2,1),(2,3,1), (3,4,1),(1,4,1)], multiedges=True)
G3.plot().show()
M3 = Matroid(G3)
display_matroid(M3, 4)

G4 = Graph([(0,1,1), (0,1,2),(1,2,1),(1,2,2), (2,3,1),(3,4,1)], multiedges=True)
G4.plot().show()
M4 = Matroid(G4)
display_matroid(M4, 4)

G5 = Graph([(0,1,1), (0,1,2),(0,1,3),(1,2,1), (2,3,1),(3,4,1)], multiedges=True)
G5.plot().show()
M5 = Matroid(G5)
display_matroid(M5, 4)


# In[7]:


# More tests
G6 = Graph([(0,1,1), (1,2,2),(0,2,1),(2,3,1), (3,4,1),(4,2,1)], multiedges=True)
G6.plot().show()
M6 = Matroid(G6)
display_matroid(M6, 4)

G7 = Graph([(0,1,1), (1,2,1),(1,2,2),(2,3,1), (3,4,1),(4,1,1)], multiedges=True)
G7.plot().show()
M7 = Matroid(G7)
display_matroid(M7, 4)


# In[8]:


def compute_fans(n):
    for i in range(0, n+1):
        if len(fan_tzfs) == i:
            z = tzf_fan(i)
            fan_tzfs.append(z)
def tzf_wheel(n):
    R.<s> = QQ['s']
    
    # cycle through some base cases ...
    if (n == 0):
        return 1
    elif (n == 1):
        return 1/(1+s)

    # prepare the fans for recursion...
    compute_fans(n-1)

    single_sum = 0
    for k in range(1, n):
        single_sum += ((-s)/(s+1))^(n-k-1) * fan_tzfs[k]
        
    # Simplified binomial sum
    non_fans = (((1-n)/(s+1)^n) * ((-s)^n-n*(-s-1)-1)) + ((n/(s+1)^n) * ((-s)^(n-1)-(n-1)*(-s-1)-1))
    
    # Note: the n = 2 triangle is a special case because the edges
    # around the wheel do not form a closed path
    if ( n > 2):
        # Brute force add cyclic matroid final case
        non_fans += ((-s)^n-1-n*(-s-1)) / ((n*s+n-1) * (s+1)^n)
    elif (n == 2):
        non_fans += 1/(1+s)

    total = n * single_sum + non_fans

    # scale the sum accordingly by the rank and size of the groundset
    if (n > 2):
        total = total * 1/(2*n*s + n)
    elif (n == 2):
        total = total / (3*s + 2)

    return total


# In[9]:


tzf_wheel(6)


# In[10]:


# Compute the toplogical zeta function for a matroid M via the
# recurrence relation summing over all proper flats
def tzf_recurrence(M):
    R.<s> = QQ['s']

    E = M.groundset()
    rk = M.full_rank()
    L = M.lattice_of_flats()
    
    # Compute base cases where the tzf is simple
    if (len(E) == 1):
        return rk * (1/(1+s))
    elif (len(E) == 0):
        return 1

    # Sum over all proper flats of L(M)
    sum_proper_flats = 0
    for F in L:
        if F == E:
            continue

        contraction = M.contract(F)
        restriction = M.delete(E - F)

        chi = x_M_reduced(contraction)
        # Recall tzf_recurrence recursively for contracted matroids
        sum_proper_flats += chi(1) * tzf_recurrence(restriction)

    return 1/(len(E) * s + rk) * sum_proper_flats


# In[14]:


tzf_recurrence(matroids.PG(0,7))


# In[111]:


def tzf_PG_calc(r, q):
    
    R.<s> = QQ['s']
    
    # Add in initial value
    tzf_PG = []
    tzf_PG.append(1/(s+1))
    
    # Iterate through as many as we need
    for i in (1..r):
        summation = prod((1-q^j) for j in (1..i))
        for k in (1..i):
            
            # Get value for reduced char poly
            poch = prod((1-q^j) for j in (1..i-k))
            summation += q_binomial(i+1, k, q) * poch * tzf_PG[k-1]
        
        # Add new value to list
        total = summation / (q_binomial(i+1, 1, q)*s+ i+1)
            
        tzf_PG.append(total)
                
    return tzf_PG


# In[112]:


tzf_PG_calc(6,2)

