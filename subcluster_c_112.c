#include "mex.h"
#include "matrix.h"
#include <math.h>

double minimos(double pry[],int fila, int colum)
{
    int i;
    double min = pry[0]; 
    for(i=0; i < fila; i++){
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
    for(i=0; i < fila; i++){
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
    for(i=0; i < fila; i++){
            if(maxim <= pry[i+colum*fila])
            {
                maxPotIndex = i;
            }
        }
    return maxPotIndex;
}

void mexFunction(int num_outputs, mxArray *outargs[], int num_inputs, const mxArray *inargs[])
{
	int i, j, rows, cols, maxPotIndex , l, findMore, numClusters,minDistSq, maxPotRatio;//,minDist; 
	double *input_matrix, *maximum, *minimum, *RAR,*sigmas,RA,*normalize,*potsVals, *dx1,*dx2, *dxSq, *d2, refPotVal, *newX, *maxPoint,*centers, minDist,maxPotVal;
     mxArray  *mat_aux[5];
    
	if ( num_inputs < 2 )
	mexErrMsgTxt("Too few input arguments!");
	else if ( num_outputs > 10 )
	mexErrMsgTxt("Too many output arguments!");

    //original matrix
    rows = mxGetM(inargs[0]); 
	cols = mxGetN(inargs[0]); 
    
    input_matrix = mxGetPr(inargs[0]);
    RAR = mxGetPr(inargs[1]);

    //Mins
    outargs[0]=mxCreateDoubleMatrix(1,cols, mxREAL);

	minimum = mxGetPr(outargs[0]);

    for(j=0; j < cols; j++){
        minimum[j] = minimos(input_matrix,rows,j);
    }

    //Maxs
    outargs[1]=mxCreateDoubleMatrix(1,cols, mxREAL);
    
    maximum = mxGetPr(outargs[1]);
    
    for(i=0; i < cols; i++){
        maximum[i] = maximos(input_matrix,rows,i);
    }
    
    //Sigmas
    RA = RAR[0];
    outargs[2]=mxCreateDoubleMatrix(1,cols, mxREAL);
    sigmas = mxGetPr(outargs[2]);
    for(j=0; j < cols; j++){
        sigmas[j] = (RA * (maximum[j] - minimum[j])) / sqrt(8.0); 
    }
    
    //Normalize the data into a unit hyperbox using the verified minX and maxX
    outargs[3]=mxCreateDoubleMatrix(rows,cols, mxREAL); 
    normalize = mxGetPr(outargs[3]);
    
    for(i=0; i < rows; i++){
        for(j=0; j < cols; j++){
            normalize[i+j*rows] = (input_matrix[i+j*rows] - minimum[j]) / (maximum[j] - minimum[j]);
        }
	}    


//    compute the initial potentials for each data point
     outargs[4] = mxCreateDoubleMatrix(1,rows, mxREAL); 
     mat_aux[0] = mxCreateDoubleMatrix(rows,cols, mxREAL); 
     d2 = mxGetPr(mat_aux[0]);
     mat_aux[1] = mxCreateDoubleMatrix(rows,cols, mxREAL); 
     dx2 = mxGetPr(mat_aux[0]);
     potsVals = mxGetPr(outargs[4]);
      for(l=0; l < rows; l++){
         potsVals[l]=0;
        for(j=0; j < rows; j++){
             d2[j]=0;
             for(i=0; i < cols; i++){
                  dx2[j+i*rows] = (normalize[l+i*rows] - normalize[j+i*rows])/RA;
                  d2[j] = d2[j] + (dx2[j+i*rows]*dx2[j+i*rows]);//^2;
             }
              potsVals[l] = potsVals[l] + exp(-4*d2[j]);     
         }
    }
    
    mat_aux[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    newX = mxGetPr(mat_aux[2]);
    
    for(i=0; i < cols; i++){
         newX[i] = referencia(potsVals,rows,cols,i);
    }
    maxPotIndex = newX[0];
    refPotVal  =  maximos(potsVals,cols,1);

//     outargs[5]=mxCreateDoubleMatrix(1,rows, mxREAL);
// 
// 	centers = mxGetPr(outargs[5]);
    
// 
//     for(i=0; i < rows; i++){
//         
//          if(i == maxPotIndex)
//             {
//                 potsVals[i]  = 0.000;
//             }
//             else
//             {
//             potsVals[i] = potsVals[i];
//             }//refPotVal;
// 
//     }
// 
// 	for(i=0; i < rows; i++){
//         
//             centers[i] = potsVals[i];
//      
//     }

// /* Start iteratively finding cluster centers and subtracting potential
//  from neighboring data points.  maxPotVal is the current highest
//  potential value and maxPotIndex is the associated data point's index. */
    
}