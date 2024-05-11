/************************************************************************************************************
* Joey Mulé, Chris Bispels																					*
* CMSC 441, Dr. Sherman																						*
* 5/10/24																									*
* Code References:																							*
* https://github.com/YYYYYW/Matrix-Multiplication															*
* https://www.geeksforgeeks.org/implementing-coppersmith-winograd-algorithm-in-java/						*
* https://github.com/Bruce-Lee-LY/matrix_multiply															*
* https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=3030fa2aecda339d593b86a260bfab9988b42df7	*
* https://www.geeksforgeeks.org/strassens-matrix-multiplication-algorithm-implementation/	                *
* Command: valgrind ./BAMspace                                                                              *
************************************************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

/*
 * BAM Algorithm
 * matA a M*K matrix
 * matB a K*N matrix
 * matC a M*N matrix
 * matC = matA * matB
 */
static void basic_matrix_multiplication(float* matA,float* matB,float* matC,const int M,const int N,const int K)
{
  for (int i = 0; i < M;i++)
    {
      for (int j = 0; j < N;j++)
	{
	  float sum = 0.0f;
	  for (int k = 0; k < K;k++)
	    {
	      sum += matA[i*K + k] * matB[k*N + j];
	    }
	  matC[i*N + j] = sum;
	}
    }
}



/*
 * Test function
 */
static void mm_test(int M, int N, int K, int rangeTop){
  unsigned seed = time(0);
  srand(seed);
  clock_t start,end;
  for(int i = 0; i < 1; i++){
    float * mA = (float*) malloc(M*K*sizeof(float));
    float * mB = (float*) malloc(K*N*sizeof(float));
    float * mD = (float*) malloc(M*N*sizeof(float));
    for(int j = 0; j < M*K; j++){
      mA[j] = rand() % rangeTop;
    }
    for(int j = 0; j < K*N; j++){
      mB[j] = rand() % rangeTop;
    }
    
    // Time stamping
    start = clock();
    basic_matrix_multiplication(mA, mB, mD, M, N, K);
    end = clock();
    double endtime = (double) (end-start)/CLOCKS_PER_SEC;
    printf("BAM%d time: %fms\n", i, endtime*1000);

    printf("\n");
    free(mA);
    free(mB);
    free(mD);
  }
}

int main(){;
  int M, N, K, rangeTop;
  M = 200;
  N = 200;
  K = 200;
  rangeTop = 100;
  printf("start\n\n");
  mm_test(M, N, K, rangeTop);
  return 0;
}
