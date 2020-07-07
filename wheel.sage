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
        summand = ((-1*s)/(s+1))^(n-k-1)*fan_tzfs[k]
        double_sum += summand
        
    non_fans = 0
    for r in range(0, n-1):
        non_fans += binomial(n, r) * (n-1-r) * (-1)^(n-1-r) * 1/(s+1)^r

    # Note: the n = 2 triangle is a special case because the edges
    # around the wheel do not form a closed path
    if ( n > 2):
        # Brute force add cyclic matroid final case
        non_fans += ((-s)^n-1-n*(-s-1)) / ((n*s+n-1) * (s+1)^n)
    elif (n == 2):
        non_fans += 1/(1+s)

    total = n * double_sum + non_fans

    # scale the sum accordingly by the rank and size of the groundset
    if (n > 2):
        total = total * 1/(2*n*s + n)
    elif (n == 2):
        total = total / (3*s + 2)

    return total

