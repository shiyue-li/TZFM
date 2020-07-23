def mobius_check(M):
    R.<s> = QQ['s']
    
    E = M.groundset()
    P = M.lattice_of_flats()
    
    total = 0
    
    for F in P:
        m = P.moebius_function(F, E)
        restrict = tzf_recurrence(M.delete(E-F))
        total += m * restrict
    
    return total
