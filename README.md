#codimensions-unstable-loci

# Extremal Values for Codimensions of Unstable Loci
Calculating extremal values for codimensions of unstable loci of SL_2 subgroups of SL_n

- Valdemar Tsanov (University of Bochum; valdemar.tsanov at rub.de)
- Yana Staneva (University of Cologne; ystaneva at math.uni-koeln.de, yanastaneva8 at gmail.com)

Date: Oct 2021

## Description:
Let _X=G/B_ be the flag variety of a semisimple complex Lie group _G_, _B_ being a fixed Borel subgroup. We are interested in the description of the GIT-classes of ample line bundles on _X_ w.r.t. a given reductive subgroup. More specifically, we are interested in extremal values of the codimensions of unstable loci. These values are obtained by a Weyl group calculation with fundamental weights, which is implemented here.

The basic step is the following: for a given pair of vectors _h_,_A_ (dominant weights), we compute the minimal length _L_ of a Weyl group element _w_ such that the scalar product _(A,wh)<0_, and we present an expression for _w_ as a product of simple reflections.

The number _L_ is interpreted suitably as the codimension of the unstable locus of the one parameter subgroup _-h_ acting on the line bundle defined by _A_.

For our purposes, we run the procedure for fixed _h_, and let _A_ run over
 - the fundamental weights
 - the sum of the fundamental weights
for every _w_. 
Our choices for _h_ are also among the above, with a few additional cases, as needed.

## Prerequisites:
Python 3.9

You need a running version of Python 3. To setup Visual Studio Code, see here: https://code.visualstudio.com/docs/python/python-tutorial

Packages needed: numpy, itertools, pandas, tabulate, tqdm

# Instructions:
1. Clone the repository
2. Open any .py script corresponding to the group type you are interested in using, e.g. E_6.py
3. Scroll down to the function called _generate_reflections(n)_ and uncomment the line which contains the number of _generators_ of the desired length _L_ you wish to test, i.e. if you want to test reflections of length _L=5_, then the _generate_reflections(n)_ function should look like:

```Python
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
```
4. You are now ready to run/debug the script.


## Comments:
The current scripts deliver results for lengths up to 10 due to hardware limitations. We strongly believe that the code can be optimized and improved for future use. Any comments and/or feedback are more than welcome.

The results and code are made available via the generous support provided by the Collaborative Research Centre/Transregio (CRC/TRR 191) on “Symplectic Structures in Geometry, Algebra and Dynamics” (http://www.mi.uni-koeln.de/CRC-TRR191/)
