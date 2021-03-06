# Assuming pwd=TZFM
# Path may need to be reconfigured for local setup
load('matroidzeta.sage')

# evaluate the first n derivatives of the tzf of matroid
# M at 0 and print the results
def get_derivatives(n, M):
	R.<s> = QQ['s']

	Z = tzf_recurrence(M)

	for i in range(0, n+1):
		print('Z^(' + str(i) + ')(0) = ' + str((-1)^i * Z(0)))
		Z = derivative(Z, s)

def compare_derivatives(n, M1, M2=None):
	R.<s> = QQ['s']

	if(M2 == None):
		E = len(M1.groundset())
		r = M1.full_rank()
		M2 = matroids.Uniform(r,E)

	Z1 = tzf_recurrence(M1)
	Z2 = tzf_recurrence(M2)

	for i in range(0, n+1):
		print('Z_M1^(' + str(i) + ')(0) = ' + str((-1)^i * Z1(0)) + 
			',\t Z_M2^(' + str(i) + ')(0) = ' + str((-1)^i * Z2(0)) + 
			'\t(' + str((-1)^i * (Z1(0)-Z2(0))) + ')')
		Z1 = derivative(Z1,s)
		Z2 = derivative(Z2,s)

def dependent_sets(M):
	E = M.groundset()
	P_E = Subsets(E, submultiset=True)

	dependent_subsets = [0] * (len(E) + 1)
	for S in P_E:
		if M.is_dependent(S):
			#print(S)
			dependent_subsets[len(S)] += 1

	for i in range(0, len(E) + 1):
		print('(' + str(dependent_subsets[i]) + ') dependent subsets of cardinality ' + str(i))

def parallel_classes(M):
	E = M.groundset()

	pc = []

	considered_edges = []

	for e in E:
		if e not in considered_edges:
			considered_edges.append(e)
			c = 1
			for f in E:
				if M.is_dependent((e,f)):
					considered_edges.append(f)
					c+= 1
			pc.append(c)

	pc.sort()
	return pc




