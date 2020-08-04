# Assuming pwd=TZFM
# Path may need to be reconfigured for local setup
load('matroidzeta.sage')

# Get the nth braid matroid
def get_braid(n):
	return Matroid(graphs.CompleteGraph(n))

# store zeta functions of braid matroids
# to dynamically compute future zeta functions
braid_tzfs = []

# compute the number of ways to choose the number of ways
# to divide n elements into a partition
def partition_choose(n, partition):
	product = 1
	consecutive = 1
	last = 0

	for part in partition:
		# choose the next part from the partition
		product = product * binomial(n, part)
		n -= part

		# discount equal partition classes
		if (part == last):
			consecutive += 1
		else:
			consecutive = 1
			last = part
		product = product / consecutive

	return product

# evaluate the reduced characteristic polynomial of
# the nth braid matroid B_n at 1
def chi_braid(n):
	return prod([(1-l) for l in range(2,n)])

# dynamically compute the topological zeta function
# of the nth braid matroid by iterating over all
# integer partitions of n
def tzf_braid(n):
	R.<s> = QQ['s']

	# cycle through some base cases ...
	if (n == 1):
		return 1
	elif (n == 2):
		return 1/(1+s)

	sum = 0
	for partition in Partitions(n):
		# sum over all proper partitions of n
		if partition == [n]:
			continue

		chi_r = chi_braid(len(partition))
		
		zetas = 1
		# dynamically compute the zeta functions in the partition
		for part in partition:
			zetas = zetas * braid_tzfs[part - 1]

		sum += partition_choose(n, partition) * chi_r * zetas

	# don't forget to multiply by the constant rank and ground set
	return 1/(binomial(n,2) * s + (n-1)) * sum

m = 20
# compute the topological zeta functions for the first m braid matroids
for i in range(1, m + 1):
	z = tzf_braid(i)
	#print(z)
	braid_tzfs.append(z)

# Computes the first m braid matroid gamma function
def braid_gamma(m):
    
    R.<s> = QQ['s']
    res = 0
    
    braids = [1, 1, -s/(s+1)]
    
    for n in (3..m):
        total = 0
        for k in (2..n):
            
            dictionary = {}
            for i in (0..n-k):
                dictionary['x{}'.format(i)] = braids[i+1]
            
            bell = bell_polynomial(n, k)
            if k != n:
                total += bell(**dictionary) * (n*(n-1)/2*s + n-k)
            else:
                total += n*(n-1)/2*s + n-k
        
        braids.append(-1/(n*(n-1)/2*s +n-1) * total)
    
    return braids
