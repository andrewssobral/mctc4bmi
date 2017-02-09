#include "mex.h"

/* Computes a part of the projection of a sparse tensor 
 * into the tangent space. See calcProjection.m for the full routine
 * THIS ROUTINE IS ONLY CALLED THROUGH CALCPROJECTION.M!
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
void calcProjection( double* result, double* result1, double* result2, double* result3, 
                     double* index, double* deriv, 
                     double* U1, double* U2, double* U3, 
                     mwSize k1, mwSize k2, mwSize k3, 
                     mwSize n1, mwSize n2, mwSize n3, 
                     mwSize nnz )
{
    mwIndex p, q, r, idx;
    mwIndex i, j, k;

    mwIndex k12;
    double temp;
    double temp2;
    
    k12 = k1*k2;

    for( idx=0; idx < nnz; ++idx ) {
        /* get the indices of the nonzero elements */
        i = index[ nDims*idx ] - 1;
        j = index[ nDims*idx + 1] - 1;
        k = index[ nDims*idx + 2] - 1;

        for( r=0; r < k3; ++r ) {
            temp = deriv[ idx ]* U3[ r + k3*k ];

            for( q=0; q < k2; ++q ) {
                temp2 = temp * U2[ q + k2*j ]; 
                for( p=0; p < k1; ++p ) {
                    result[ AT3(p,q,r) ] += temp2 * U1[ p + k1*i ];
                }
         
                result1[ i + n1*q + n1*k2*r ] += temp2;
            }

            for( p=0; p < k1; ++p ) {
                result2[ p + j*k1 + k1*n2*r ] += temp * U1[ p + k1*i ];
            }
        }

        for( q=0; q < k2; ++q ) {
            for( p=0; p < k1; ++p ) {
                result3[ AT3(p,q,k) ] += deriv[ idx ] * U1[ p + k1*i ] * U2[ q + k2*j ];
            }
        }
    }
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double* deriv;
    double* U1;
    double* U2;
    double* U3;

    double* result;               /* output tensor */
    double* result1;              /* output tensor */
    double* result2;              /* output tensor */
    double* result3;              /* output tensor */

    double* index;
    mwSize nnz;
    const mwSize* dims;
    mwSize k1, k2, k3;
    mwSize n1, n2, n3;

    /* GET THE INDEX ARRAY */
    /* ------------------- */

    index = (double*) mxGetPr( prhs[0] );
    nnz = mxGetN( prhs[0] );

    if ( mxGetM(prhs[0]) != nDims )
        mexErrMsgIdAndTxt( "arrayProduct:Dimensions", 
                           "Error in first input. This function currently only works for 3D tensors.");
        

    /* GET THE VALUES OF EUCLID. DERIV*/
    /* ------------------- */

    deriv = mxGetPr( prhs[1] );
    if ( mxGetM( prhs[1] ) != nnz || 
         mxGetN( prhs[1] ) != 1 )
        mexErrMsgIdAndTxt( "arrayProduct:Dimensions", 
                           "Error in second input. Must be a (number of nonzeros) column vector.");
    
    /* GET THE U FACTORS       */ 
    /* ----------------------- */
    
    U1 = mxGetPr( prhs[2] );
    dims = mxGetDimensions( prhs[2] );
    k1 = dims[0];
    n1 = dims[1];

    U2 = mxGetPr( prhs[3] );
    dims = mxGetDimensions( prhs[3] );
    k2 = dims[0];
    n2 = dims[1];

    U3 = mxGetPr( prhs[4] );
    dims = mxGetDimensions( prhs[4] );
    k3 = dims[0];
    n3 = dims[1];

    /* create the output vector (vectorized result tensor) */
    plhs[0] = mxCreateDoubleMatrix( k1*k2*k3, 1, mxREAL );
    plhs[1] = mxCreateDoubleMatrix( n1*k2*k3, 1, mxREAL );
    plhs[2] = mxCreateDoubleMatrix( n2*k1*k3, 1, mxREAL );
    plhs[3] = mxCreateDoubleMatrix( n3*k1*k2, 1, mxREAL );
    
    /* get a pointer to the real data in the output matrix */
    result = mxGetPr( plhs[0] );
    result1 = mxGetPr( plhs[1] );
    result2 = mxGetPr( plhs[2] );
    result3 = mxGetPr( plhs[3] );

    /* call the computational routine */
    calcProjection( result, result1, result2, result3, index, deriv, U1, U2, U3, k1, k2, k3, n1, n2, n3, nnz );
}
