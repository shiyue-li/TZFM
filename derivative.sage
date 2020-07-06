# Assuming pwd=TZFM
# Path may need to be reconfigured for local setup
load('matroidzeta.sage')

def get_derivatives(M, n):
	R.<s> = QQ['s']

	Z = tzf_recurrence(M)

	for i in range(0, n+1):
		print('Z^(' + str(i) + ')(0) = ' + str((-1)^i * Z(0)))
		Z = derivative(Z, s)

