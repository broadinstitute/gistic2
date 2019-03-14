/* Mex/C inplementation of cosort 2012-09-05 Steven Schumacher DFCI */
/*
  A Matlab function for supporting the SegArray class:

      [group val index] = cosort(v1,...,vN)

  where inputs v1,...,vN are a series of already sorted double vectors (typically indices
  of SegArray breakpoints
*/
#include <string.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) 
{
  if (nlhs>3)
    mexErrMsgTxt("maximum number of outputs for cosort is 3.");

  /* get size of input lists */
  int comb_size = 0;
  int *vsizes = (int *) mxCalloc(nrhs,sizeof(int));
  double **vptrs = mxCalloc(nrhs,sizeof(double*));
  int *vidx = mxCalloc(nrhs,sizeof(int));
  int i;
  for (i=0; i<nrhs; i++) {
    if( mxIsDouble(prhs[i]) ) {
      vsizes[i] = mxGetNumberOfElements(prhs[i]);
      comb_size += vsizes[i];
      vptrs[i] = mxGetPr(prhs[i]);
      vidx[i] = 0;
    }
    else {
      mexErrMsgTxt("non-double index passed to cosort()");
    }
  }

  /* create output structures */
  mxArray* group = mxCreateDoubleMatrix(comb_size,1,mxREAL);
  double* pgroup = mxGetPr(group);
  plhs[0] = group;

  mxArray* sorted = NULL; /* Matlab handle for optional "sorted" value output */
  double* psorted = NULL; /* C pointer for "sorted" output */

  if( nlhs >= 2 ) {
    sorted = mxCreateDoubleMatrix(comb_size,1,mxREAL);
    plhs[1] = sorted;
    psorted = mxGetPr(sorted);
  }
  mxArray* index = NULL; /* Matlab handle for optional "index" output */
  double* pindex = NULL; /* C pointer for "index" output */
  if( nlhs >= 3 ) {
    index = mxCreateDoubleMatrix(comb_size,1,mxREAL);
    plhs[2] = index;
    pindex = mxGetPr(index);
  }
  if( comb_size > 0 ) {
    int k;
    for( k=0; k<comb_size; k++ ) {
      /* find minimum from all the sorted lists */
      double minv = 0; /* minimum value */
      int mini = 0;    /* index of minimum value */
      int i;
      bool gotval = 0;
      for( i=0; i<nrhs; i++ ) {
	if( vidx[i] < vsizes[i] ) { 
	  if( gotval ) {
	    double val = *vptrs[i];
	    if( val < minv ) {
	      mini = i;
	      minv = val;
	    }
	  }
	  else {
	    mini = i;
	    minv = *vptrs[i];
	    gotval = 1;
	  }
	}
      }
      /* update input pointers */
      vptrs[mini]++;
      vidx[mini]++;

      /* update output(s) element */
      pgroup[k] = (double) (mini+1);
      if( nlhs >= 2 )
	psorted[k] = (double) minv;
      if( nlhs >= 3 )
	pindex[k] = (double) vidx[mini];
    }
  }
  /* free allocated blocks */
  mxFree(vsizes);
  mxFree(vptrs);
  mxFree(vidx);

}


