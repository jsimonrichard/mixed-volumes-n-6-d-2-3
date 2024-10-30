import itertools 

conf1 = 2* [(0,0,1),(1,0,1),(0,1,1)]
conf2 = [(a,b,1) for a in [0,1] for b in [0,1]] + [(1,0,0),(0,1,0)]


def permuted_configurations(configs):
    """
        all permutations of the point configurations 
        from the list in the input 
    """
    set_perm = lambda c : set(itertools.permutations(c))
    result = set()
    for c in configs: 
        result.update(set_perm(c))
    return result 

def prod_map_abs_gr(c):
    M = matrix(c)
    assert M.nrows() == 6 and M.ncols() == 3
    return tuple([abs(M[list(I),:].det()*M[list(J),:].det()) for I,J in SetPartitions(range(6),[3,3])])

rays = set([prod_map_abs_gr(c) for c in permuted_configurations([conf1,conf2])]) 

    
C_6_3 = Polyhedron(rays = rays)
print("Indexing:")
for i,I in enumerate(SetPartitions(range(6),[3,3])):
    print("{}: {} {}".format(i+1,list(I[0]),list(I[1])))
    
print()

ineqs = [ineq[1:] for ineq in C_6_3.inequalities_list() ]
print(matrix(ineqs))

