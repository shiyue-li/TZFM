# Assuming pwd=TZFM
# Path may need to be reconfigured for local setup
load('matroidzeta.sage')

# Get the nth thagomizer graph, K_1,1,n
def get_thag_graph(n):
	# Construct edge set
	top_edges = [(0,i) for i in range(2, n+2)]
	bottom_edges = [(1,i) for i in range(2, n+2)]
	edges = top_edges + bottom_edges + [(0,1)]

	g = Graph(edges) 
	return g

# Get the nth thagomizer matroid
def get_thag_matroid(n):
	g = get_thag_graph(n)
	T_n = Matroid(g)
	return T_n

# Helper functions to computer the closed form of the
# TZF for thagomizer matroids 

# Compute the product (2k + 1)s + (k+1) over specified range
# Appearing in closed form of TZF for thagomizer matroids
def thag_product(low_bound, up_bound):
	R.<s> = QQ['s']
	product = 1
	for k in range(low_bound, up_bound + 1):
		product = product * ((2*k + 1)*s + k+ 1)

	return product

# Compute the sum n!/k! * (s+1)^(n-k) (1-s)^k * sum (2j + 1)s + (j+1)
# Appearing in closed form of TZF for thagomizer matroids
def thag_sum(n):
	R.<s> = QQ['s']
	sum = 0

	for k in range(1, n + 1):
		summand = (factorial(n))/(factorial(k))
		summand = summand * ((s+1)^(n-k)) * ((1-s)^k)
		summand = summand * thag_product(1, k-1)
		sum = sum + summand

	return sum

# Compute the closed form of the TZF for the nth thagomzier matroid
def tzf_thag(n):
	R.<s> = QQ['s']

	numerator = factorial(n) * ((s+1)^(n-1)) + thag_sum(n)
	denominator = ((s+1)^n) * thag_product(1,n)

	tzf = numerator / denominator
	return tzf

# Check that the closed form of the TZF for thagomizer matroids 
# Gives the correct function for T_n
def check_thag_tzf(n):
	closed_form = tzf_thag(n)

	T_n = get_thag_matroid(n)
	definition = tzf(T_n)

	difference = closed_form - definition

	# Not actually sure if this works ...
	# Sage seems to have trouble verifying if
	# different rational expressions are equal

	print(difference)

check_thag_tzf(6)