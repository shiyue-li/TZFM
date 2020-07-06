#!/usr/bin/env python
# coding: utf-8

# In[46]:


# Create some fan graphs
list_of_fans = []
n=5
for m in range(0, n+1):
    edges = [(0,i) for i in range(1, m+1)] + [(i, i+1) for i in range(1, m)]
    list_of_fans.append(Matroid(Graph(edges)))

# Get lattices
lattices_of_fans = [F.lattice_of_flats() for F in list_of_fans]    
    
def get_fan(n):
    # Make fan graph F(1,n) using edges
    edges = [(0,i) for i in range(1, n+1)] + [(i, i+1) for i in range(1, n)]
    
    # Make graph
    g = Graph(edges)    
    F = Matroid(g)
    F_L = F.lattice_of_flats()
    
    # Get level sets of lattice
    level_sets = F_L.level_sets()
    sizes = [len(level) for level in level_sets]

    print('F(1,{0}):'.format(n), sizes)

# 8 is as feasible as it gets
for i in range(1, 6):
    get_fan(i)


# In[24]:


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


# In[3]:


# MAX CODE FOR TZF
# given a lattice L and a flag F in L containing at least the min and max
# of L, xmf(L,F) is the product of reduced characteristic polynomials of
# intevals in F, and xmf1(L,F) finds the value of this polynomial when you set
# q = 1.

def xmf(L,F):
    l = len(F)
    R.<q> = QQ['q']
    numerator = product([posetify(L.interval(F[i],F[i+1])).characteristic_polynomial() for i in [0..(l-2)]])
    denominator = (q - 1)^(l-1)
    return (numerator/denominator)

def xmf1(L,F):
    l = len(F)
    R.<q> = QQ['q']
    numerator = product([posetify(L.interval(F[i],F[i+1])).characteristic_polynomial() for i in [0..(l-2)]])
    denominator = (q - 1)^(l-1)
    return (numerator/denominator)(1)

# given a lattice L and a flag F in L containing at least the min and
# max of L, trf(L,F) finds the rational function (in the variable s)
# that you need when building the topological zeta function.
def trf(L,F):
    R.<s> = QQ['s']
    rk = L.rank_function()
    return product([1/((len(f))*s + rk(f)) for f in F if not f == L.bottom()])
    
# Utility: For a set S of subsets, form its associated poset w.r.t. containment
def posetify(S):
    return Poset((S,[[a,b] for a in S for b in S if a.issubset(b)]))

def tzf(M):
    L = M.lattice_of_flats()
    flags = [c for c in L.chains() if len(c) > 0 and c[0] == L.bottom() and c[-1] == L.top()]
    return sum([xmf1(L,f)*trf(L,f) for f in flags])


# In[29]:


# Matroid relax test
M = matroids.named_matroids.NonVamos()
print(M)
MR = relax_circuit_hyperplane(M)

print(tzf(M))
print(tzf(MR))


# In[34]:


# Get some TZFs
for n in range(1,6):

    # Fan graph
    edges = [(0,i) for i in range(1, n+1)] + [(i, i+1) for i in range(1, n)]
    g_f = Graph(edges)    
    print('Fan({0}): '.format(n), tzf(Matroid(g_f)))
    
    # Complete graph
    # g_k = graphs.CompleteGraph(n+1)    
    # print('Complete({0}): '.format(n+1), tzf(Matroid(g_k)))

