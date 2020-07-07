# Assuming pwd=TZFM
# Path may need to be reconfigured for local setup
load('matroidzeta.sage')

# Get the nth fan graph: one vertex adjoined 
# to a path with n vertices
def get_fan_graph(n):
	# Construct edge set
	spokes = [(0,i) for i in range(1, n+1)]
	path = [(i, i+1) for i in range(1,n)]
	edges = spokes + path

	g = Graph(edges) 
	return g

# Get the nth fan matroid
def get_fan_matroid(n):
	g = get_fan_graph(n)
	F_n = Matroid(g)
	return F_n

# store zeta functions of fan matroids
# to dynamically compute future zeta functions
fan_tzfs = []

# dynamically compute the topological zeta function
# of the nth fan matroid by summing over smaller fans
def tzf_fan(n):
	R.<s> = QQ['s']

	# cycle through some base cases ...
	if (n == 0):
		return 1
	elif (n == 1):
		return 1/(1+s)

	# summing over pairs of smaller fans
	tzf = sum(fan_tzfs[k] * fan_tzfs[n - 1 - k] for k in (0..n-1))
	tzf = tzf * 1/((2*n-1)*s+n)

	# contribution of previous fan...
	tzf += ((2*n-3)*s + (n-1))/((2*n-1)*s+n) * (-s/(s+1)) * fan_tzfs[n-1]

	return tzf

m = 20
# compute the topological zeta functions for the first m fan matroids
for i in range(0, m + 1):
	z = tzf_fan(i)
	#print(z)
	fan_tzfs.append(z)