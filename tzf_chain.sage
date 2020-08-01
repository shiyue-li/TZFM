# Calls a function on the components and multiplies them together
def componentize(func):
    def wrapper(*args, **kwargs):
        if not args or not isinstance(args[0], sage.matroids.matroid.Matroid):
            return 'Componentize needs a function that takes a matroid as the first argument!'
        
        # Get components and call function on them
        components = args[0].components()
        E = args[0].groundset()
        prod = 1
        for c in components:
            prod *= func(Matroid(args[0] \ (E - c)), *args[1:], **kwargs)
        return prod
    return wrapper

# Compute the toplogical zeta function for a matroid M, using chains
@componentize
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
