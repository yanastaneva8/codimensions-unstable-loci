
# coding: utf-8

# In[ ]:


"""Extremal Values for Codimensions of Unstable Loci
Authors: Valdemar Tsanov, Yana Staneva
Date: August 24, 2020

Let X=G/B be the flag variety of a semisimple complex Lie group G, B being a fixed Borel subgroup. 
We are interested in the description of the GIT-classes of ample line bundles on X w.r.t. a given semisimple subgroup

This script computes the extremal values for codimensions of unstable loci.
SL_2-subgroups of SL_n
Type A_n

Python 3.7
"""


# Import necessary python packages, if unavailable, use command - pip3 install numpy
import numpy as np
# Import necessary python packages, if unavailable, use command - pip3 install itertools
import itertools
from numpy.linalg import inv


# Initialize a list to store the computations for the main loop
comp_list = []
# Defines a class to store all computations
class Computations: 
    def __init__(self, binary, h, length, weyl, scalar_product):
        self.binary = binary
        self.h = h
        self.length = length
        self.weyl = weyl
        self.scalar_product = scalar_product

    def __str__(self):
        return "binary: {}, h: {}, length: {}, weyl: {}, scalar_product: {} ".format(self.binary, self.h, self.length, self.weyl, self.scalar_product)
    
    def __getitem__(self, item):
        return getattr(self, item)


# Defines a function to generate the fundamental weights b
def binary_a(n):
    # Initialize an empty nx1 array 
    b = np.zeros(n)
    # Initialize an identity nxn array  
    binary = np.zeros(shape=(n,n))
    # Iterate over all values up to n
    for i in range(0,n):
        # Insert a 1 in the i-th entry of the nx1 vector b
        np.put(b, i, 1)
        # Insert the resulting b vector in the initialized array to store
        binary[i]=b   
    return binary


# Defines a function to generate the vectors h
def h_vector(n):
    ## Name of M?
    # Generate the nxn matrix M, s.t. h = Mv, v is a vector with entries 1 & 0
    M = np.triu(2*np.ones(n))
    M[:, -1] = np.ones(n)
    # print('M = %s' %M)
    # Call function for fundamental weights b to obtain all possible vectors v s.t. h = Mv
    b = binary_a(n)
    # Initialize empty list to store the resulting vectors h
    h_array=[]
    # Iterate over all possible fundamental weights b
    for i in range(len(binary_a(n))):
        # Compute the dot product h = Mv for v=b_i
        h = np.dot(M, np.array([b[i]]).T)
        # Store resulting vector h in the list h_array
        h_array.append(np.asarray(h))
        # print('h_array = %s' % np.ravel(h)   ) 
    return h_array


# Defines a function to generate the permutation in matrix form w.r.t. the simple reflections of the group
def matrix_form(elt):
    # print(elt)
    # Assign the symbolic expression of the simple reflections, e.g. r_1,...,r_n, r_1*r_2,...,r_1*r_n, ...
    # and create a list which stores the symbolic expression generated above  
    elts = list(elt)
    # print(elts)
    # Store only the indeces of the generated reflections, i.e. if w=r_1*r_2*r_3, only store 123
    reflections = elts[1::3]
    # print(reflections)
    # Initialize empty nxn identity matrix
    matrixform = np.eye(n)
    # Loop to iterate over all generated reflection expressions
    for elt in reflections:
        # Assign the index of the current reflection to the variable a
        a = int(elt)
        # Initialize empty nxn identity matrix to generate the basis elements
        mat = np.eye(n)
        # The n-th simple reflection element with -1 in the (n,n) position
        if a == n:
            mat[a-1, a-1] = -1
            # Compute the matrix form for the n-th simple reflection
            matrixform = np.dot(matrixform, mat)
        # The rest j=1,...,n-1 simple reflections for type A_n, the block ((0,1),(1,0)) shifted along the diagonal to begin at j-th row    
        else:
            mat[a-1, a] = 1
            mat[a, a - 1] = 1
            mat[a,a] = 0
            mat[a-1, a-1] = 0
            # Compute the matrix form for the j-th simple reflection
            matrixform = np.dot(matrixform, mat)
    return matrixform


# Initialize two lists, one for the generating symbolic reflections (r_1,r_2,...,r_j,...,r_n), 
# one for the tuples of reflections, i.e. w = r_1*r_2*...*r_j.
# To reduce complexity, comment in or out the lines which generate the n-tuples of reflections needed,
# i.e. if you want to test length l=3, then leave uncommented lines 132-133
generators = []
all_possible_roots = []
def generate_reflections(n):
    for i in range(1, n+1):
        generators.append("r"+str(i))
    # The line below generates expressions with length 1
    all_possible_roots.extend(generators)   
    # The line below generates expressions with length 2
    for subset in itertools.product(generators, generators):
        all_possible_roots.append('*'.join(subset)) 
    # The line below generates expressions with length 3
    for subset in itertools.product(generators, generators, generators):
        all_possible_roots.append('*'.join(subset))
    # The line below generates expressions with length 4
    for subset in itertools.product(generators, generators, generators, generators):
        all_possible_roots.append('*'.join(subset))
    # The line below generates expressions with length 5
#     for subset in itertools.product(generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # The line below generates expressions with length 6
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # The line below generates expressions with length 7
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # The line below generates expressions with length 8
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # The line below generates expressions with length 9
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # The line below generates expressions with length 10
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    return(all_possible_roots)


# Input here the dimension
n = 3
# print('dimension n = %s' %n)
# Call function to generate the vectors h
h_store = h_vector(n)
# print('array with h vectors, h_store = %s' %h_store)
# Call function to generate fundamental weights b
binary = binary_a(n)
# print('fundamental weights, b_i = %s' %binary)
# Call function to generate reflections w
generators = generate_reflections(n)
# print('reflection, symbolic_w = %s' %generators)
# Initialize list to store the reflections which yield a nonzero scalar product
nonzero_reflections = [] 
# Initialize list to store the reflections in matrix form which yield a nonzero scalar product
nonzero_reflections_matrix_form = []
# Initialize a list to store the resulting scalar products
scalar = []
# Main loop to go over all permutations generated
for i in generators:
    # print('reflection symbolic expression, w = %s' %i)
    # Convert the generated reflection into matrix form
    a = matrix_form(i)
    # print('reflection in matrix form, w = %s' %a)
    # Store the current reflection
    nonzero_reflections.append(i)
    # Store the current reflection in matrix form
    nonzero_reflections_matrix_form.append(a)
    # Nested loop to iterate all fundamental weights b
    for j in range(len(binary)):
        # print('fundamental weight, b = %s' %binary[j])
        # Apply the generated permutation on the fundamental weight in the current loop
        weyl_acting_b = np.dot(a, binary[j])
        # print('reflection applied to binary vector, wb = %s' %weyl_acting_b)
        # Nested loop to go over all generated h vectors
        for h in h_store:
            # Compute the scalar product of the applied permutation onto the fundamental weight, wb, and the vector h
            scalar_prod = np.dot(weyl_acting_b, h)
            # Conditional statement for storing the scalar product
            if scalar_prod<0:
                # print('generator expression = %s' %i)
                # print('matrix form = %s' %matrix_form(i))
                # print('binary = %s' %binary[j])
                # print('weyl acting on binary = %s' % weyl_acting_b)
                # print(' h_vector = %s' %h)
                # print('scalar = %s' %scalar_prod)
                # Store the current scalar product if it satisfies the condition
                scalar.append(scalar_prod)
                # Call the class to store the computations
                temp_permutation = Computations(binary[j], np.ravel(h), i.count("r"), i, scalar_prod)
                # Append the current computation if satisfying the condition
                comp_list.append(temp_permutation)
# printing tests           
# print(scalar)
# print(len(scalar))
# print('r1 = %s' % matrix_form('r1'))
# print('r2 = %s' % matrix_form('r2'))
# print('r3 = %s' % matrix_form('r3'))
# print('r4 = %s' % matrix_form('r4'))
        
for item in comp_list:
    print(item)



""" Script to convert the resulting computations into a LaTeX table
"""


# Import necessary python packages, if unavailable, use command - pip3 install pandas
import pandas as pd
# Import necessary python packages, if unavailable, use command - pip3 install tabulate
import tabulate
from tabulate import tabulate


# Defines a function to print output from class Computations in ready to paste LaTeX table format
def latex_with_lines(df, *args, **kwargs):
    kwargs['column_format'] = '|'.join([''] + ['l'] * df.index.nlevels
                                            + ['r'] * df.shape[1] + [''])
    kwargs['index'] = False
    res = df.to_latex(*args, **kwargs)
    return res.replace('\\\\\n', '\\\\ \\hline\n')


# Call function to print output in ready to paste in LaTeX table format
df = pd.DataFrame({'Ah': x.h.astype(int), 'Binary': x.binary.astype(int), 'Length': x.length, 'Scalar': x.scalar_product, 'CWeyl': x.weyl} for x in comp_list)
df_to_print = latex_with_lines(df)
print(df_to_print)
#print([x.weyl_matrix for x in comp_list])
#tabulate([x.weyl_matrix for x in comp_list], tablefmt="latex")
#print(tabulate((x.weyl_matrix.astype(int) for x in comp_list), tablefmt="latex"))

