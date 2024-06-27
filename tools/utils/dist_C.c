#include "mex.h"

/*DIST_C Distances between vectors.
 *	
 *	dist_C(X,Y)
 *	  X - L1xR matrix of rows vectors.
 *	  Y - RxL2 matrix of column vectors.
 *	
 *  if wrong number of arguments or matrix sizes do not match
 *
 *  Returns Error
 *	
 *  Else
 *
 *  Returns L1xL2 matrix of vectors distances.
 *
 *  For optimal results X should have many more rows than columns and Y
 *  less columns than X rows
 *  This can be done in MatLab by manipulating input and output:
 *  dist_C(X,Y) == dist_C(Y',X')'
 *
 *  This is a MEX-file for MATLAB.
 *  By Robert Lau, RWTH-Student 2014.
 */

#define X           prhs[0] //should be L1xR matrix
#define Y           prhs[1] //should be RxL2 matrix
#define Z           plhs[0] //should become L1xL2 matrix

void mexFunction(int nlhs, mxArray *plhs[], /* Output variable */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
  size_t r, l1, l2, i, j, k;
  double *p_Z, *p_X, *p_Y;

  if(nrhs != 2)
    mexErrMsgTxt("Wrong number of input arguments");
  if(nlhs != 1)
    mexErrMsgTxt("Too many output arguments");

  r = mxGetN(X);
  if(r != mxGetM(Y))
    mexErrMsgTxt("Matrix sizes do not match");
  l1 = mxGetM(X);
  l2 = mxGetN(Y);

  Z = mxCreateDoubleMatrix(l1, l2, mxREAL);
  p_Z = mxGetPr(Z);
  p_X = mxGetPr(X);
  p_Y = mxGetPr(Y);
  if(r == 1) {
    for(i = 0; i < l1; i++) {
      for(j = 0; j < l2; j++) {
        p_Z[j*l2 + i] = fabs(p_X[i] - p_Y[j]);
      }
    }
  } else {
    for(i = 0; i < l1; i++) {
      for(j = 0; j < r; j++) {
        for(k = 0; k < l2; k++) {
          p_Z[k*l1 + i] += (p_X[i + j*l1] - p_Y[r*k + j]) * (p_X[i + j*l1] - p_Y[r*k + j]);
          if(j == r-1) p_Z[k*l1 + i] = sqrt(p_Z[k*l1 + i]);
        }
      }
    }
  }
}