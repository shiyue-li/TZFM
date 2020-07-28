load('matroidzeta.sage')

#The Gamma function
def Gamma(M,k=0):
    R.<s> = QQ['s']
    
    L = M.lattice_of_flats()
    E = M.groundset()
    # Sum over all proper flats of L(M)
    
    sum_proper_flats = 0
    for F in L:
        
        restriction = M.delete(E - F)
        Z = tzf_recurrence(restriction)
        
        contraction = M.contract(F)
        chi = x_M(contraction)
        chi_0 = chi(q=0)
        
        if type(chi) is Integer:
            sum_proper_flats += chi *Z
        else:
            sum_proper_flats += chi_0 *Z
        print(F,'---',M.rank(F))

    return sum_proper_flats

# I had problems getting the above implementation fo gamma to work for me...
def gamma(M, k=0):
    R.<s> = QQ['s']
    L = M.lattice_of_flats()
    E = M.groundset()
     
    sum_flats = 0
    for F in L:
        restrict = M.delete(E-F)
        Z = tzf_recurrence(restrict)
        Z_der = derivative(Z, s, k)
        mob = L.moebius_function(F, E)
        sum_flats += mob * Z_der
    return sum_flats

def gamma_recurrence(M):
    R.<s> = QQ['s']

    E = M.groundset()
    rk = lambda S : M.rank(S)
    L = M.lattice_of_flats()
    
    # Compute base cases where the tzf is simple
    if (rk(E) == 1):
        return -len(E)*s/(len(E)*s+1)
    elif (rk(E) == 0):
        return 1

    # Sum over all proper flats of L(M)
    sum_flats = 0
    for F in L:
        if F == E:
            continue

        contraction = M.contract(F)
        restriction = M.delete(E - F)

        # Recall tzf_recurrence recursively for contracted matroids
        sum_flats += (len(E)*s + rk(F)) * gamma_recurrence(restriction)

    return -1/(len(E) * s + rk(E)) * sum_flats

# some work in progress ideas...
def gamma_f(M, f=tzf_recurrence):
    R.<s> = QQ['s']
    L = M.lattice_of_flats()
    E = M.groundset()
     
    sum_flats = 0
    for F in L:
        restrict = M.delete(E-F)
        Z = f(restrict)
        mob = L.moebius_function(F, E)
        sum_flats += mob * Z
    return sum_flats


#A test to see what recurrence relation Gamma satisfies
def rcurr_gamma(M):  
    R.<s> = QQ['s']
    E = M.groundset()
    L = M.lattice_of_flats()
    sum_ = 0
    for F in L:
        if len(F) == 0:
            #sum_+= -s/(s+1)
            continue
        if F == E:
            continue
        restriction = M.delete(E-F)
        contraction = M.contract(F)
        sum_+= 1/(len(F)*s+M.rank(F))*(Gamma(restriction))
        #print(F, '----', sum_)
    
    return sum_

