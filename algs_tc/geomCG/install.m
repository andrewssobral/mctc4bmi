%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

disp('                ___   ___      |    geomCG: Tensor completion by Riemannian Optimization')
disp('   _  _  _     / __|/ __|      |    ')
disp('  / \/_\/ \|\/| (__ |(_/\      |    Daniel Kressner')
disp('  \_|\__\_/|  |\___|\___/      |    Michael Steinlechner')
disp('  \_|                          |    Bart Vandereycken')
disp('     ')
disp('geomCG is licensed under a BSD 2-clause license, see LICENSE.txt')
disp('For an explanation of the algorithm, we refer to')
disp('   ')
disp('      D. Kressner, M. Steinlechner, B. Vandereycken:')
disp('      <a href="http://sma.epfl.ch/~anchpcommon/publications/tensorcompletion.pdf">Low-Rank Tensor Completion by Riemannian Optimization</a>.')
disp('      MATHICSE Technical Report 20.2013, June 2013. Submitted to BIT Numerical Mathematics.')

disp('   ')
                                   
disp('Compiling mex files...')
mex calcFunction_mex.c
mex calcGradient_mex.c
mex calcInitial_mex.c
mex calcProjection_mex.c
mex getValsAtIndex_mex.c

addpath( cd )
addpath( [cd, filesep, 'examples'] )
disp('Make sure that the Tensor Toolbox is also properly installed and loaded into the current path')
disp('You can download the Tensor Toolbox from <a href="http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html">Sandia Labs</a>. GeomCG was tested with version 2.5 only!')

disp('Finished. Try out the example code simpleGeomCG.m')









