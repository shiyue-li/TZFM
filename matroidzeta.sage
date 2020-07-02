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

# given a lattice L and a flag F in L containing at least the min and
# max of L, mrf(L,F) finds the rational function (in the variable s)
# that you need when building the motivic zeta function.
def mrf(L,F):
    R.<q,T> = QQ['q,T']
    rk = L.rank_function()
    return product([(T^(len(f))*(q-1))/(q^(rk(f)) - T^(len(f))) for f in F if not f == L.bottom()])
    
# Utility: For a set S of subsets, form its associated poset w.r.t. containment
def posetify(S):
    return Poset((S,[[a,b] for a in S for b in S if a.issubset(b)]))

def tzf(M):
    L = M.lattice_of_flats()
    flags = [c for c in L.chains() if len(c) > 0 and c[0] == L.bottom() and c[-1] == L.top()]
    return sum([xmf1(L,f)*trf(L,f) for f in flags])

def mzf(M):
    R.<q,T> = QQ['q,T']
    L = M.lattice_of_flats()
    n = len(M.groundset())
    flags = [c for c in L.chains() if len(c) > 0 and c[0] == L.bottom() and c[-1] == L.top()]
    return (1/T^n)*sum([xmf(L,f)*mrf(L,f) for f in flags])
