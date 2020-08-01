# Copute the toplogical zeta function for a matroid M, using chains
def tzf_chain(M):
    
    R.<s> = ZZ['s']
    res = 0
    
    # Get lattice of flats and chains
    L = M.lattice_of_flats()
    chains = [x for x in list(L.chains()) if x and x[0] == frozenset()]

    # Calculate ranks
    ranks = {}
    rk = L.rank_function()
    for ele in L.list():
        ranks[ele] = L.rank(ele)
            
    # Run over each chain
    for c in chains:
        
        prod = 1
        
        # Do multiplications
        for i in (1..len(c)-1):
            prod *= (-len(c[i])*s - ranks[c[i-1]])
        
        # Do divisions (for some reason it is faster to do them separately)
        for i in (1..len(c)-1):
            prod /= (len(c[i])*s + ranks[c[i]])
        
        res += prod
    
    return res
