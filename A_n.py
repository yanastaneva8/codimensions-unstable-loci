
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
"""


import numpy as np
import itertools
from numpy.linalg import inv

# initialize a list to store the computations for the main loop
comp_list = []
# defines a class to store all computations
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


# defines a function to generate the fundamental weights b
def binary_a(n):
    b = np.zeros(n)
    binary = np.ones(shape=(n,n))
    for i in range(n):
        np.put(b, i, 1)
        binary[i]=b
    return binary


# defines a function to generate the vectors h
def h_vector(n):
    M = np.triu(2*np.ones(n))
    M[:, -1] = np.ones(n)
    3# print('M = %s' %M)
    b = binary_a(n)
    h_array=[]
    for i in range(len(binary_a(n))):
        h = np.dot(M, np.array([b[i]]).T)
        h_array.append(np.asarray(h))
        # print('h_array = %s' % np.ravel(h)   ) 
    return h_array

# defines a function to generate the permutation in matrix form w.r.t. the basis elements of the group
def matrix_form(elt):
    elts = list(elt)
    reflections = elts[1::3]
    matrixform = np.eye(n)
    for elt in reflections:
        a = int(elt)
        mat = np.eye(n)
        if a == n:
            mat[a-1, a-1] = -1
            matrixform = np.dot(matrixform, mat)
        else:
            #mat[a - 2, a - 2] = 1
            mat[a-1, a] = 1
            mat[a, a - 1] = 1
            mat[a,a] = 0
            mat[a-1, a-1] = 0
            matrixform = np.dot(matrixform, mat)
    return matrixform


# initialize two lists, one for the generating expressions of the permutations, one for the the list of expressions
# when needed, comment in or out the lines which generate the n-tuples of permutations
generators = []
all_possible_roots = []
def generate_reflections(n):
    for i in range(1, n+1):
        generators.append("r"+str(i))
    # the line below generates expressions with length 1
    all_possible_roots.extend(generators)   
    # the line below generates expressions with length 2
    for subset in itertools.product(generators, generators):
        all_possible_roots.append('*'.join(subset)) 
    # the line below generates expressions with length 3
    for subset in itertools.product(generators, generators, generators):
        all_possible_roots.append('*'.join(subset))
    # the line below generates expressions with length 4
    for subset in itertools.product(generators, generators, generators, generators):
        all_possible_roots.append('*'.join(subset))
    # the line below generates expressions with length 5
#     for subset in itertools.product(generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # the line below generates expressions with length 6
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # the line below generates expressions with length 7
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # the line below generates expressions with length 8
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # the line below generates expressions with length 9
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    # the line below generates expressions with length 10
#     for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators, generators, generators):
#         all_possible_roots.append('*'.join(subset))
    return(all_possible_roots)

# input here the dimension
n = 3
# call function to generate the vectors h
h_store = h_vector(n)
#print('h_store = %s' %h_store)
# call function to generate fundamental weights b
binary = binary_a(n)
#print('b = %s' %binary)
# call function to generate reflections r
generators = generate_reflections(n)
#print('reflection = %s' %generators)
# initialize list to store the reflections which yield a nonzero scalar product
nonzero_reflections = [] 
# initialize list to store the reflections in matrix form which yield a nonzero scalar product
nonzero_reflections_matrix_form = []
# initialize a list to store the resulting scalar products
scalar = []
# Main loop to go over all permutations generated
for i in generators:
    #print('refl = %s' %i)
    # Convert the generated reflection into matrix form
    a = matrix_form(i)
    #print('reflmatrix = %s' %a)
    # Store the current reflection
    nonzero_reflections.append(i)
    # Store the current reflection in matrix form
    nonzero_reflections_matrix_form.append(a)
    # Nested loop to go over the list of fundamental weights
    for j in range(len(binary)):
        #print('b = %s' %binary[j])
        # Apply the generated permutation on the fundamental weight in the current loop
        weyl_acting_b = np.dot(a, binary[j])
        #print('weyl = %s' %weyl_acting_b)
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

    
# [1 1 , ... 0] implement 


# In[ ]:


""" Script to convert the resulting computations into a LaTeX table
"""


import pandas as pd
import tabulate
from tabulate import tabulate


def latex_with_lines(df, *args, **kwargs):
    kwargs['column_format'] = '|'.join([''] + ['l'] * df.index.nlevels
                                            + ['r'] * df.shape[1] + [''])
    kwargs['index'] = False
    res = df.to_latex(*args, **kwargs)
    return res.replace('\\\\\n', '\\\\ \\hline\n')

df = pd.DataFrame({'Ah': x.h.astype(int), 'Binary': x.binary.astype(int), 'Length': x.length, 'Scalar': x.scalar_product, 'CWeyl': x.weyl} for x in comp_list)
df_to_print = latex_with_lines(df)
print(df_to_print)
#print([x.weyl_matrix for x in comp_list])
#tabulate([x.weyl_matrix for x in comp_list], tablefmt="latex")
#print(tabulate((x.weyl_matrix.astype(int) for x in comp_list), tablefmt="latex"))

