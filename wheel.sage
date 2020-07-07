# Assuming pwd=TZFM
# Path may need to be reconfigured for local setup
load('matroidzeta.sage')

# Get the nth wheel graph: here we adopt the convention
# that the nth wheel graph has n "spokes"
def get_wheel_graph(n):
	# Construct edge set
	spokes = [(0,i) for i in range(1, n+1)]
	path = [(i, i+1) for i in range(1,n)]

	if (n > 1):
		path += [(1,n)]

	edges = spokes + path

	g = Graph(edges) 
	return g

# Get the nth wheel matroid
def get_wheel_matroid(n):
	g = get_wheel_graph(n)
	W_n = Matroid(g)
	return W_n

load('fan.sage')

def compute_fans(n):
	for i in range(0, n+1):
		if len(fan_tzfs) == i:
			z = tzf_fan(i)
			fan_tzfs.append(z)

def tzf_wheel(n):
	R.<s> = QQ['s']

	# cycle through some base cases ...
	if (n == 0):
		return 1
	elif (n == 1):
		return 1/(1+s)

	# prepare the fans for recursion...
	compute_fans(n-1)

	double_sum = 0
	for k in range(1, n):
		for r in range(0, n-k):
			double_sum += binomial(n-1-k, r) * (-1)^(n-r-k) * 1/(s+1)^r * fan_tzfs[k]

	non_fans = 0
	for r in range(0, n-1):
		non_fans += binomial(n, r) * (n-1-r) * (-1)^(n-1-r) * 1/(s+1)^r

	total = n * double_sum + non_fans

	return 1/(2*n*s + n) * total

