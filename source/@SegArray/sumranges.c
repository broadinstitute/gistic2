/* Mex/C inplementation of sumranges 2012-09-05 Steven Schumacher DFCI */
/*
  A Matlab function for supporting the SegArray class:

      values = sumranges(B,E,V,size)

*/
#include <string.h>
#include "mex.h"
/*#include "matrix.h"
 */

void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) 
{
  /* input validation */
  if ( nlhs > 1 )
    mexErrMsgTxt("sumranges() returns at most one output");
  if ( nrhs != 4 )
    mexErrMsgTxt("sumranges() requires four arguments");
  if ( 1 != mxGetNumberOfElements(prhs[3]) )
    mexErrMsgTxt("sumranges() argument 4 must be a scalar");

  /* get output size (argument 4) */
  double *poutput_size = mxGetPr(prhs[3]);
  int Nout = (int)*poutput_size;
  /* allocate and initialize array for output */
  mxArray* outvals = mxCreateDoubleMatrix(Nout,1,mxREAL);
  double* poutvals = mxGetPr(outvals);
  plhs[0] = outvals;
  int i;
  for( i=0; i<Nout; i++) {
    poutvals[i] = 0.0;
  }

  /* process vector arguments */
  int Nin = mxGetNumberOfElements(prhs[0]); 
  double* pStarts = mxGetPr(prhs[0]);
  if( Nin != mxGetNumberOfElements(prhs[1]) || Nin != mxGetNumberOfElements(prhs[2]) )
    mexErrMsgTxt("all of the first 3 arguments of sumranges() must have the same length");
  double* pEnds = mxGetPr(prhs[1]);
  double* pValues = mxGetPr(prhs[2]);

  for( i=0; i<Nin; i++) {
    int start = ((int)pStarts[i])-1;
    int end = ((int)pEnds[i])-1;
    if( start >= Nout || end >= Nout )
      mexErrMsgTxt("sumranges index exceeds size ouf output");
    int j;
    for( j=start; j<=end; j++) {
      poutvals[j] += pValues[i];
    }
  }
}


