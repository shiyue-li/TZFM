
# truncate a matroid M to rank r
def truncate(M, r):
    ind_sets = []
    for i in range(0, r+1):
        ind_sets += M.independent_r_sets(i)
    trunc = Matroid(M.groundset(), independent_sets=ind_sets)
    return trunc

# here is some code I wrote to check my work for the
# girth - zeta function conjecture

# computes the elemntary symmetric functions on x_i = i
# this is a cool dynamic algorithm, I like it
# perhaps it could be used in some kind of induction proof...?
def elem_symmetric(k):
    x = [i for i in range(1,k)]
    a = [1] + [0] * (k-1)
    for n in range(0, k-1):
         for j in reversed(range(1,k)):
             a[j] = a[j] + a[j-1] * x[n]
    return a

# computes stirling numbers of the second kind
def stirling(n, k):
    sum = 0
    for j in range(0, k+1):
        sum += (-1)^(k-j) * binomial(k, j) * j^n
    return 1 / factorial(k) * sum

# the malevolent sum from the conjecture...
# 1 < j < k
def check_sum(k, j):
    sum = 0
    a = elem_symmetric(k-1)
    for l in range(j, k):
        sum += stirling(l,j) * ((k-j) * a[k-l] - (k-1)*j*a[k-l-1])
    return sum + (k-j)*a[0]*stirling(k,j)