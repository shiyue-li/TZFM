
# truncate a matroid M to rank r
def truncate(M, r):
     ind_sets = []
     for i in range(0, r+1):
         ind_sets += M.independent_r_sets(i)
     trunc = Matroid(M.groundset(), independent_sets=ind_sets)
     return trunc