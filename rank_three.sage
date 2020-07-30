load('gamma.sage')
import itertools

def findsubsets(s, n): 
    return list(itertools.combinations(s, n)) 

def convert(s): 
    str1 = ""    
    return(str1.join(s))

# get the size four circuits of a rank 3 simple matroid
# i.e. find all four-element subsets of E that do not 
# contain a given size three circuit
def get_four_circuits(ground_set, three_circuits):
    four_circuits = []
    
    four_subsets = findsubsets(ground_set, 4)
    for S in four_subsets:
        contains_circuit = False
        for C in three_circuits:
        	# check if C \subseteq S
            contains_circuit = contains_circuit or all(x in S for x in C)
        if not contains_circuit:
        	# this isn't a circuit!
            four_circuits.append(convert(S))
    return four_circuits

# find a rank 3 (simple) matroid by inputting the ground set
# and size three circuits -> i.e. given a ground set and a 
# list of lines, find me the rank 3 matroid
def get_simple_3_matroid(ground_set, three_circuits):
    four_circuits = get_four_circuits(ground_set, three_circuits)
    circuits = three_circuits + four_circuits
    return Matroid(groundset=ground_set, circuits=circuits)

# These are the two matroids in Max's paper that are not isomorphic
# but have the same zeta function. This calculation
# shows that they actually have the same gamma function!
E = 'abcdefg'
N1 = get_simple_3_matroid(E, ['abe','adg','bcd','efg'])
N2 = get_simple_3_matroid(E, ['abe','adg','ecd','efg'])

# Note: the labeling is the same as in Max's paper
# N1 has the "horizontal" fourth circuit, while N2 has the "diagonal" fourth circuit

# They are not isomorphic matroids
N1.is_isomorphic(N2)

# But they have the same zeta function...
Z_N1 = tzf_recurrence(N1)
Z_N2 = tzf_recurrence(N2)

# And the same gamma function...
T_N1 = gamma_recurrence(N1)
T_N2 = gamma_recurrence(N2)

# here are their lattices (sans labeling) ...
N1.lattice_of_flats().plot(label_elements=False, title="L(N1)").show()
N2.lattice_of_flats().plot(label_elements=False, title="L(N2)").show()