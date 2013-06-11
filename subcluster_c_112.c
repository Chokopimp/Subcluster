#include "mex.h"
#include "matrix.h"
#include <math.h>

double minimos(double pry[],int fila, int colum)
{
  int i;
  double min = pry[0]; 
  for(i=0; i < fila; i++)
  {
    if(min > pry[i+colum*fila])
    {
      min = pry[i+colum*fila];
    }
  }
  return min ;
}

double maximos(double pry[],int fila, int colum)
{
  int i;
  double max = pry[0+colum*fila]; 
  for(i=0; i < fila; i++)
  {
    if(max < pry[i+colum*fila])
    {
      max = pry[i+colum*fila];
    }
  }
  return max;  
}

double referencia (double pry[],int fila, int column, int colum)
{
  int i;
  double maxim, maxPotIndex;
  maxim = maximos(pry,fila,colum);
  for(i=0; i < fila; i++)
  {
    if(maxim <= pry[i+colum*fila])
    {
      maxPotIndex = i;
    }
  }
  return maxPotIndex;
}

void mexFunction(int num_outputs, mxArray *outargs[], int num_inputs, const mxArray *inargs[])
{
	int i, j, rows, cols, arar,maxPotIndex , l, findMore, numClusters, maxPotRatio;//,minDist; 
  double *input_matrix, *maximum, *minimum, *RAR,*sigmas,RA,*normalize,*potVals , *dx1,*dx2, dxSq,minDistSq, *d2, refPotVal, *newX, *maxPoint,*centers, minDist,maxPotVal;
  mxArray  *mat_aux[6];
    
	if ( num_inputs < 2 )
    mexErrMsgTxt("Too few input arguments!");
	else if ( num_outputs > 10 )
    mexErrMsgTxt("Too many output arguments!");

    //original matrix
  rows = (int)mxGetM(inargs[0]); 
	cols = (int)mxGetN(inargs[0]); 
    
  input_matrix = mxGetPr(inargs[0]);
  RAR = mxGetPr(inargs[1]);

    //Mins
  outargs[0]=mxCreateDoubleMatrix(1,cols, mxREAL);

	minimum = mxGetPr(outargs[0]);

  for(j=0; j < cols; j++)
  {
    minimum[j] = minimos(input_matrix,rows,j);
  }

    //Maxs
  outargs[1]=mxCreateDoubleMatrix(1,cols, mxREAL);
    
  maximum = mxGetPr(outargs[1]);
    
  for(i=0; i < cols; i++)
  {
    maximum[i] = maximos(input_matrix,rows,i);
  }
    
    //Sigmas
  RA = RAR[0];
  outargs[2]=mxCreateDoubleMatrix(1,cols, mxREAL);
  sigmas = mxGetPr(outargs[2]);

  for(j=0; j < cols; j++)
  {
    sigmas[j] = (RA * (maximum[j] - minimum[j])) / sqrt(8.0); 
  }
    
  //Normalize the data into a unit hyperbox using the verified minX and maxX
  mat_aux[0]=mxCreateDoubleMatrix(rows,cols, mxREAL); 
  normalize = mxGetPr(mat_aux[0]);
    
  for(i=0; i < rows; i++)
  {
    for(j=0; j < cols; j++)
    {
      normalize[i+j*rows] = (input_matrix[i+j*rows] - minimum[j]) / (maximum[j] - minimum[j]);
    }
	}    

//    compute the initial potentials for each data point
  mat_aux[1] = mxCreateDoubleMatrix(1,rows, mxREAL); 
  potVals = mxGetPr(mat_aux[1]);
  mat_aux[2] = mxCreateDoubleMatrix(rows,cols, mxREAL); 
  d2 = mxGetPr(mat_aux[2]);
  mat_aux[3] = mxCreateDoubleMatrix(rows,cols, mxREAL); 
  dx2 = mxGetPr(mat_aux[3]);
 

  for(l=0; l < rows; l++)
  {
    potVals[l]=0;
    for(j=0; j < rows; j++)
    {
      d2[j]=0;
      for(i=0; i < cols; i++)
      {
        dx2[j+i*rows] = (normalize[l+i*rows] - normalize[j+i*rows])/RA;
       d2[j] = d2[j] + (dx2[j+i*rows]*dx2[j+i*rows]);//^2;
      }
      potVals[l] = potVals[l] + exp(-4*d2[j]);     
    }
  }

  mat_aux[4] = mxCreateDoubleMatrix(1,1, mxREAL);
  newX = mxGetPr(mat_aux[4]);
    
  for(i=0; i < cols; i++)
  {
    newX[i] = referencia(potVals,rows,cols,i);
  }
  maxPotIndex = (int)newX[0];
  refPotVal  =  maximos(potVals,cols,1);
  maxPotVal = refPotVal;

 /* Start iteratively finding cluster centers and subtracting potential
  from neighboring data points.  maxPotVal is the current highest
  potential value and maxPotIndex is the associated data point's index. */
  outargs[3]=mxCreateDoubleMatrix(rows, cols, mxREAL);
  centers = mxGetPr(outargs[3]);
  numClusters = 0;
  findMore = 1;
  mat_aux[5] = mxCreateDoubleMatrix(1,cols, mxREAL);
  maxPoint = mxGetPr(mat_aux[5]);

  while(findMore && maxPotVal)
  {
    findMore = 0;
    
    for(i=0; i<cols; i++)
    {
      maxPoint[i] = normalize[maxPotIndex+i*rows];
    }

    maxPotRatio = maxPotVal/refPotVal;
    
    if(maxPotRatio > 0.5)
    {
      findMore =1; //the new peak value is significant, accept it
    }
    
    else if(maxPotRatio > 0.15) /*accept this data point only if it achieves a good balance between having a reasonable potential and being far from all existing cluster centers*/
    {
      minDistSq = -1;

      for(j=0; j<numClusters; j++)
      {
        dxSq = 0.0;
        for(i=0; i<cols ;i++)
        {
          dx1[i]=(maxPoint[i] - centers[j+i*rows])/RA;
          dxSq = dxSq + (dx1[i]* dx1[i]);// dx1[i]^2
        }

        if(minDistSq < 0 || dxSq < minDistSq)
        {
          minDistSq = dxSq;
        }
      }

      minDist = sqrt(minDistSq);

      if((maxPotRatio + minDist) > 1)
      {
        findMore = 1;
      }

      else
      {
        findMore = 2;
      }
    }
    
    if(findMore == 1)
    {
      

      for(i=0; i<cols; i++)
      {
        centers[(numClusters)+i*rows] = maxPoint[i];
      }

      for(j = 0; j < rows; j++) //subtract potential from data points near the new cluster center.
      {
        d2[j] = 0.0;
       
        for(i = 0; i < cols; i++ )
        {
          dx2[j+i*rows] = (maxPoint[i]-normalize[j+i*rows])/(1.25*RA);
          d2[j] = d2[j] + (dx2[j+i*rows]*dx2[j+i*rows]);
        }
        potVals[j] = potVals[j] - maxPotVal*exp(-4*d2[j]);//deduct(pp);
       
        if(potVals[j]<0)
        {
          potVals[j] = 0;
        }
      }

      for(i=0; i < cols; i++)
      {
        newX[i] = referencia(potVals,rows,cols,i);
      }

      maxPotIndex = newX[0]; 
      refPotVal  =  maximos(potVals,cols,1);
      maxPotVal = refPotVal; //[maxPotVal,maxPotIndex] = maximos(potVals); //find the data point with the highest remaining potential
      numClusters = numClusters + 1;
    }
    
    else if(findMore == 2)
    {
      potVals[maxPotIndex] = 0;
      
      for(i=0; i < cols; i++)
      {
        newX[i] = referencia(potVals,rows,cols,i);
      }

      maxPotIndex = newX[0]; 
      refPotVal  =  maximos(potVals,cols,1);
      maxPotVal = refPotVal; //[maxPotVal,maxPotIndex] = maximos(potVals);
    }
  }

  for(j=0; j<numClusters; j++)
  {
    for(i=0; i<cols; i++)
    {
      centers[j+i*numClusters] = (centers[j+i*numClusters] * (maximum[i] - minimum[i])) + minimum[i];
    }
  }
}