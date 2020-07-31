def Gamma_flag(M):
    R.<s> = QQ['s']
    
    L = M.lattice_of_flats()
    flags = [c for c in L.chains() if len(c) > 0 and c[0] == L.bottom() and c[-1] == L.top()]
    for F in flags:
        l = len(F)
        prod =1
        prod*= (-1)^(l-1)*product([(len(F[i])*s+M.rank(F[i-1]))/(len(F[i])*s+M.rank(F[i])) for i in [1..(l-1)]])
        #print(prod,'----', F)
