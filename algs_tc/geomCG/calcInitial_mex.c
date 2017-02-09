#include "mex.h"

/* Computes the initial guess for the line search.
 * THIS ROUTINE IS ONLY CALLED THROUGH CALCINITIAL.M!
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
void calcFunction( double* result, double* index, double* R, 
                   double* S, double* U1, double* U2, double* U3, 
                   double* Y_tilde, double* U1_tilde, double* U2_tilde, double* U3_tilde, 
                   mwSize k1, mwSize k2, mwSize k3, 
                   mwSize nnz )
{
    mwIndex p, q, r, idx;
    mwIndex i, j, k;

    mwIndex k12;

    double A1, A2,
           B1, B2,
           C1, C2, 
           D2, sum, sum2;
    
    k12 = k1*k2;

/* CALCULATION OF:
 *    res = Y_tilde x U1 x U2 x U3 + S x U1_tilde x U2 x U3 + S x U1 x U2_tilde x U3 + ...
 *          ( A                  ) + ( B                  ) + ( C                  ) + ( D ...
 */
    sum2 = 0.;

    for( idx=0; idx < nnz; ++idx ) {
        /* get the indices of the nonzero elements */
        i = index[ nDims*idx ] - 1;
        j = index[ nDims*idx + 1] - 1;
        k = index[ nDims*idx + 2] - 1;

        A1 = 0.;
        A2 = 0.;

        B1 = 0.;
        B2 = 0.;

        C1 = 0.;
        C2 = 0.;

        D2 = 0.;

        sum = 0.;

        for( r=0; r < k3; ++r ) {
            for( q=0; q < k2; ++q ) {
                for( p=0; p < k1; ++p ) {
                    
                    A1 += Y_tilde[ AT3(p,q,r) ] * U1[ p + k1*i ];
                    C1 += S[ AT3(p,q,r) ]       * U1[ p + k1*i ]; 
                    /* D1 is same as C1, thus not needed */

                    B1 += S[ AT3(p,q,r) ]       * U1_tilde[ p + k1*i ];
                }

                A2 += A1 * U2[ q + k2*j ];
                B2 += B1 * U2[ q + k2*j ];
                D2 += C1 * U2[ q + k2*j ];

                C2 += C1 * U2_tilde[ q + k2*j ];

                A1 = 0.;
                B1 = 0.;
                C1 = 0.;
            }

            sum += ( A2 + B2 + C2 ) * U3[ r + k3*k ]
                              + D2 * U3_tilde[ r + k3*k ];

            A2 = 0.;
            B2 = 0.;
            C2 = 0.;
            D2 = 0.;
        }
        sum2 += sum * sum;
        *result += sum * R[idx];
    }
    *result /= sum2;
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double* R;
    double* S;
    double* U1;
    double* U2;
    double* U3;

    double* Y_tilde;
    double* U1_tilde;
    double* U2_tilde;
    double* U3_tilde;

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
        

    /* GET THE VALUES OF EUCLID DERIV */
    /* ------------------- */

    R = mxGetPr( prhs[1] );
    if ( mxGetM( prhs[1] ) != nnz || 
         mxGetN( prhs[1] ) != 1 )
        mexErrMsgIdAndTxt( "arrayProduct:Dimensions", 
                           "Error in second input. Must be a (number of nonzeros) column vector.");

    /* GET THE CORE TENSOR     */
    /* ----------------------- */

    S = mxGetPr( prhs[2] );

    dims = mxGetDimensions( prhs[2] );
    k1 = dims[0];
    k2 = dims[1];

    if (nDims == mxGetNumberOfDimensions( prhs[2] ))
        k3 = dims[2];
    else
        k3 = 1;

    
    /* GET THE U FACTORS       */ 
    /* ----------------------- */
    
    U1 = mxGetPr( prhs[3] );
    U2 = mxGetPr( prhs[4] );
    U3 = mxGetPr( prhs[5] );



    Y_tilde = mxGetPr( prhs[6] );

    /*if (nDims != mxGetNumberOfDimensions( prhs[6] ))
        mexErrMsgIdAndTxt( "arrayProduct:Dimensions", 
                           "Error in sixth input. This function currently only works for 3D tensors.");*/
    U1_tilde = mxGetPr( prhs[7] );
    U2_tilde = mxGetPr( prhs[8] );
    U3_tilde = mxGetPr( prhs[9] );

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );

    /* get a pointer to the real data in the output */
    result = mxGetPr( plhs[0] );

    /* call the computational routine */
    calcFunction( result, index, R, S, U1, U2, U3, Y_tilde, U1_tilde, U2_tilde, U3_tilde, k1, k2, k3, nnz );
}
