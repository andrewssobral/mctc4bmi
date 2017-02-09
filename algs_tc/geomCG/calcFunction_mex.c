#include "mex.h"

/* Computes the value of the objective function
 *
 *    || X_Omega - A_Omega || ^2
 *
 * between a sparse matrix A and
 * a Tucker tensor X.
 * THIS ROUTINE IS ONLY CALLED THROUGH CALCFUNCTION.M!
 * DO NOT CALL SEPARATELY!
 *
 * GeomCG Tensor Completion. Copyright 2013 by
 * Michael Steinlechner
 * Questions and contact: michael.steinlechner@epfl.ch
 * BSD 2-clause license, see LICENSE.txt
 * 
*/

/* HELPER FUNCTIONS FOR LINEAR ARRAY ACCESS */
#define AT2(i,j) i+n*j
#define AT3(i,j,k) i+k1*j+k12*k

/* GLOBAL VARIABLE DEFINITION FOR THE NUMBER OF DIMENSIONS */
const mwSize nDims = 3;


/* The computational routine */
void calcFunction( double* result, double* index, double* A, 
                   double* G, double* U1, double* U2, double* U3, 
                   mwSize k1, mwSize k2, mwSize k3, 
                   mwSize nnz )
{
    mwIndex p, q, r, idx;
    mwIndex i, j, k;

    mwIndex k12;
    double temp;
    double temp2;
    double temp3;
    
    k12 = k1*k2;

    for( idx=0; idx < nnz; ++idx ) {
        /* get the indices of the nonzero elements */
        i = index[ nDims*idx ] - 1;
        j = index[ nDims*idx + 1] - 1;
        k = index[ nDims*idx + 2] - 1;

        temp = 0.;
        temp2 = 0.;
        temp3 = 0.;
        for( r=0; r < k3; ++r ) {
            for( q=0; q < k2; ++q ) {
                for( p=0; p < k1; ++p ) {
                    
                    temp += G[ AT3(p,q,r) ] * U1[ p + k1*i ] ;
                
                }
                temp2 += temp*U2[ q + k2*j ];
                temp = 0.;

            }
            temp3 += temp2*U3[ r + k3*k ];
            temp2 = 0.;
        }

        temp3 -= A[idx];
        *result += temp3 * temp3;
    }
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double* A;
    double* G;
    double* U1;
    double* U2;
    double* U3;

    double* result;              /* output matrix */

    double* index;
    mwSize nnz;
    const mwSize* dims;
    mwSize k1, k2, k3;

    /* GET THE INDEX ARRAY */
    /* ------------------- */

    index = (double*) mxGetPr( prhs[0] );
    nnz = mxGetN( prhs[0] );

    if ( mxGetM(prhs[0]) != nDims )
        mexErrMsgIdAndTxt( "arrayProduct:Dimensions", 
                           "Error in first input. This function currently only works for 3D tensors.");
        

    /* GET THE VALUES OF A */
    /* ------------------- */

    A = mxGetPr( prhs[1] );
    if ( mxGetM( prhs[1] ) != nnz || 
         mxGetN( prhs[1] ) != 1 )
        mexErrMsgIdAndTxt( "arrayProduct:Dimensions", 
                           "Error in second input. Must be a (number of nonzeros) column vector.");

    /* GET THE CORE TENSOR     */
    /* ----------------------- */

    G = mxGetPr( prhs[2] );

    dims = mxGetDimensions( prhs[2] );
    k1 = dims[0];
    k2 = dims[1];
    
    if (nDims == mxGetNumberOfDimensions( prhs[2] ))
        k3 = dims[2];
    else
        k3 = 1;
        /*mexErrMsgIdAndTxt( "arrayProduct:Dimensions", 
                           "Error in third input. This function currently only works for 3D tensors.");*/

    
    
    

    
    /* GET THE U FACTORS       */ 
    /* ----------------------- */
    
    U1 = mxGetPr( prhs[3] );

    U2 = mxGetPr( prhs[4] );

    U3 = mxGetPr( prhs[5] );

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );

    /* get a pointer to the real data in the output */
    result = mxGetPr( plhs[0] );

    /* call the computational routine */
    calcFunction( result, index, A, G, U1, U2, U3, k1, k2, k3, nnz );
}
