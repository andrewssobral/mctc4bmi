                ___   ___   |  geomCG: 
   _  _  _     / __|/ __|   |  Tensor completion by Riemannian Optimization 
  / \/_\/ \|\/| (__ |(_/\   |  
  \_|\__\_/|  |\___|\___/   |  Michael Steinlechner
  \_|                       |  Ecole Polytechnique Federale de Lausanne


Description
-----------

geomCG is an implementation of a Riemannian tensor completion algorithm [1]. It
performs a geometric non-linear conjugate gradient optimization scheme on the
manifold of tensors of fixed multi-linear rank. 
For a detailed algorithmic description, we refer to the technical report.

It makes use of the Tensor Toolbox 2.5 [2,3], which has to be installed and added to your
Matlab path first. 

geomCG is licensed under a 2-clause BSD license, see LICENSE.txt. Feel free to use
it for your own projects. If you do so, please cite the corresponding technical
report [1].


Installation
------------

After having downloaded and installed the Tensor Toolbox 2.5, navigate to the geomCG
folder within Matlab and run install.m. A working mex compiler is NECESSARY for the
optimized C routines handling the sparse tensors.

geomCG was sucessfully tested on 
    Linux:   Matlab R2009b and R2012a
    Mac OS:  Matlab R2012a and R2013a
    Windows: Matlab R2012b

Earlier versions may not work.

Example code for geomCG is found in the subdirectory examples.


Contact
-------

For questions and suggestions contact

Michael Steinlechner
MATHICSE, Ecole Polytechnique Federale de Lausanne
michael.steinlechner@epfl.ch


References
----------

[1]   D. Kressner, M. Steinlechner, B. Vandereycken:
      Low-Rank Tensor Completion by Riemannian Optimization.
      MATHICSE Technical Report 20.2013, June 2013. Submitted to BIT Numerical Mathematics.

[2]   B. W. Bader, T. G. Kolda and others. MATLAB Tensor Toolbox Version 2.5, January 2012.
      Available from http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html

[3]   B. W. Bader and T. G. Kolda. Algorithm 862: MATLAB tensor classes for fast
      algorithm prototyping, ACM Transactions on Mathematical Software
      32(4):635-653, December 2006.
