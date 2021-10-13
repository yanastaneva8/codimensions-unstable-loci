# E_7

"""Extremal Values for Codimensions of Unstable Loci
Authors: Valdemar Tsanov, Yana Staneva
Date: Oct 1, 2021

Let X=G/B be the flag variety of a semisimple complex Lie group G, B being a fixed Borel subgroup. 
We are interested in the description of the GIT-classes of ample line bundles on X w.r.t. a given reductive subgroup. 
More specifically, we are interested in extremal values of the codimensions of unstable loci.

This script computes the extremal values for codimensions of unstable loci for one-paramater subgroups.
The one-parameter subgroups we have implemented here are:
- fundamental coweights
- the sum of the fundamental weights (or one half times the sum of the positive roots)

Type E_7

Python 3.9
"""

# Import packages, if unavailable, use command: pip3 install numpy
import numpy as np
import itertools
import pandas as pd
import tabulate
from tabulate import tabulate
# Import package for loop progress bar in console
import tqdm
from tqdm import tqdm

# Class to store computations
comp_list = []
class Computations:
    def __init__(self, fund, h, length, weyl, scalar_product):
        self.fund = fund
        self.h = h
        self.length = length
        self.weyl = weyl
        self.scalar_product = scalar_product

    def __str__(self):
        return "h: {}, fund: {}, length: {}, weyl: {}, scalar_product: {} ".format(self.h, self.fund, self.length, [int(s) for s in self.weyl if s.isdigit()], self.scalar_product)
    
    def __getitem__(self, item):
        return getattr(self, item)

# Function to store the fundamental weights in an nxn array
def fundamental(n):
    fundamental = np.zeros(shape=(n-1,n))
    fundamental[0] = np.array([0,0,0,0,0,0,-1,1]) 
    fundamental[1] = np.array([1,1,1,1,1,1,-1,1])
    fundamental[2] = np.array([-1,1,1,1,1,1,-3,3])
    fundamental[3] = np.array([0,0,1,1,1,1,-2,2])
    fundamental[4] = np.array([0,0,0,2,2,2,-3,3])
    fundamental[5] = np.array([0,0,0,0,1,1,-1,1])
    fundamental[6] = np.array([0,0,0,0,0,2,-1,1])
    return fundamental

# Function to compute the vectors h
def h_vector(n):
    funds = fundamental(n)
    h = np.vstack([funds, [0,1,2,3,4,5,-17/2,-17/2]])
    # print('h_array = %s' % h   ) 
    return h

# Function to represent the simple reflections in matrix form
def matrix_form(elt):
    elts = list(elt)
    reflections = elts[1::3]
    matrixform = np.eye(n)
    for elt in reflections:
        a = int(elt)
        mat = np.eye(n)
        if a == 1:
            mat = -1*np.ones((n,n))
            mat[a-1, ::-1] = 1
            mat[:, a-1] = 1
            mat[n-1, :] = 1
            mat[:, n-1] = 1
            mat[a-1, n-1] = -1
            mat[n-1, a-1] = -1
            np.fill_diagonal(mat, 3)
            matrixform = np.dot(0.25*matrixform, mat)
        elif a == 2:
            mat[a-1, a-2] = -1
            mat[a-2, a-1] = -1
            mat[a-2, a-2] = 0
            mat[a-1, a-1] = 0
            matrixform = np.dot(matrixform, mat)
        elif a == 8:
            mat[a-2, a-2] = 0
            mat[a-3, a-2] = 1
            mat[a-2, a-3] = 1
            mat[a-3, a-3] = 0
            matrixform = np.dot(matrixform, mat)
        else:
            mat[a-2, a-2] = 0
            mat[a-3, a-2] = 1
            mat[a-2, a-3] = 1
            mat[a-3, a-3] = 0
            matrixform = np.dot(matrixform, mat)
    return matrixform

# Initialize two lists, one for the generating symbolic reflections (r_1,r_2,...,r_j,...,r_n), 
# one for the tuples of reflections, i.e. w = r_1*r_2*...*r_j.
# To reduce complexity, comment in or out the lines which generate the n-tuples of reflections needed
generators = []
all_possible_roots = []
def generate_reflections(n):
    for i in range(1, n+1):
        generators.append("r"+str(i))
    # all_possible_roots.extend(generators)    
    # for subset in itertools.product(generators, generators):
    # for subset in itertools.product(generators, generators, generators):
    # for subset in itertools.product(generators, generators, generators, generators):
    for subset in itertools.product(generators, generators, generators, generators, generators):
    # for subset in itertools.product(generators, generators, generators, generators, generators, generators):
    # for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators):
    # for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators):
    # for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators, generators):
    # for subset in itertools.product(generators, generators, generators, generators, generators, generators, generators, generators, generators, generators):
        all_possible_roots.append('*'.join(subset))
    return(all_possible_roots)

# Function to print results in LaTeX table
def latex_with_lines(df, *args, **kwargs):
    kwargs['column_format'] = '|'.join([''] + ['l'] * df.index.nlevels
                                            + ['r'] * df.shape[1] + [''])
    kwargs['index'] = False
    res = df.to_latex(*args, **kwargs)
    return res.replace('\\\\\n', '\\\\ \\hline\n')


# # Main function to test all fundamental weights
if __name__ == "__main__":
    n = 8
    fund = fundamental(n)
    generators = generate_reflections(n-1)
    h_vectors = h_vector(n)
    nonzero_reflections = [] 
    nonzero_reflections_matrix_form = []
    scalar = []
    for i in tqdm(generators):
        a = matrix_form(i)
        nonzero_reflections.append(i)
        nonzero_reflections_matrix_form.append(a)
        for j in range(len(fund)):
            weyl_acting_fund = np.dot(a, fund[j])
            for h in fund:
                scalar_prod = np.dot(weyl_acting_fund, h)
                if scalar_prod<0:
                    scalar.append(scalar_prod)
                    temp_permutation = Computations(fund[j], np.ravel(h), i.count("r"), i, scalar_prod)
                    comp_list.append(temp_permutation)  
    # for item in comp_list:
    #     print(item)
    df = pd.DataFrame({'h': x.h.astype(int), '\varpi': x.fund.astype(int), 'Length': x.length, 'Scalar': x.scalar_product, 'Weyl': [int(s) for s in x.weyl if s.isdigit()]} for x in comp_list)
    df_to_print = latex_with_lines(df)
    print(df_to_print)
    pass

# # Main function to test specific fundamental weight
# if __name__ == "__main__":
#     n = 8
#     fund = np.array([0,0,0,0,0,2,-1,1]) 
#     generators = generate_reflections(n-1)
#     # h_vectors = h_vector(n)
#     nonzero_reflections = [] 
#     nonzero_reflections_matrix_form = []
#     scalar = []
#     for i in tqdm(generators):
#         a = matrix_form(i)
#         nonzero_reflections.append(i)
#         nonzero_reflections_matrix_form.append(a)
#         weyl_acting_fund = np.dot(a, fund)
#         scalar_prod = np.dot(weyl_acting_fund, fund)
#         if scalar_prod<0:
#             scalar.append(scalar_prod)
#             temp_permutation = Computations(fund, fund, i.count("r"), i, scalar_prod)
#             comp_list.append(temp_permutation)  
#     # for item in comp_list:
#     #     print(item)
#     df = pd.DataFrame({'h': x.fund.astype(int), '\varpi': x.fund.astype(int), 'Length': x.length, 'Scalar': x.scalar_product, 'Weyl': [int(s) for s in x.weyl if s.isdigit()]} for x in comp_list)
#     df_to_print = latex_with_lines(df)
#     print(df_to_print)
#     pass
