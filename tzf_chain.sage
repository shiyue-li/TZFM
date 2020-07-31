# Copute the toplogical zeta function for a matroid M, using chain recurrence
def tzf_chain(M):
    
    res = 0
    
    # Get lattice of flats and chains
    L = M.lattice_of_flats()
    chains = [x for x in list(L.chains()) if x and x[0] == frozenset()]

    # Calculate coefficient from rank function
    R.<s> = QQ['s']
    rk = L.rank_function()
    
    # Run over each chain
    for c in chains:
        prod = 1
        for i in (1..len(c)-1):
            prod *= -(len(c[i])*s + rk(c[i-1])) / (len(c[i])*s + rk(c[i]))
        res += prod
    
    return res
