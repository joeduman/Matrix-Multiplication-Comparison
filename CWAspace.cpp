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
* Command: valgrind ./CWAspace                                                                              *
************************************************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

/*
 * Used in Coppersmith-Winograd Algorithm
 * strideA is the col num of matA, initial value is K
 * strideB is the col num of matB, initial value is N
 * strideC is the col num of matC, initial value is N
 */
static void mm_basic(float* matA, float* matB, float* matC, const int M, const int N, const int K,
		     const int strideA, const int strideB, const int strideC) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      float sum = 0.0f;
      for (int k = 0; k < K; k++) {
	sum += matA[i * strideA + k] * matB[k * strideB + j];
      }
      matC[i * strideC + j] = sum;
    }
  }
}

/*
 * matA M*K
 * matB K*N
 * matC M*N
 * matC = matA * matB
 * S1 = A21 + A22     T1 = B12 - B11
 * S2 = S1 - A11      T2 = B22 - T1
 * S3 = A11 - A21     T3 = B22 - B12
 * S4 = A12 - S2      T4 = T2 - B21
 * M1 = A11 * B11     U1 = M1 + M2
 * M2 = A12 * B21     U2 = M1 + M6
 * M3 = S4 * B22      U3 = U2 + M7
 * M4 = A22 * T4      U4 = U2 + M5
 * M5 = S1 * T1       U5 = U4 + M3
 * M6 = S2 * T2       U6 = U3 - U4
 * M7 = S3 * T3       U7 = U3 + M5
 * C11 = U1
 * C12 = U5
 * C21 = U6
 * C22 = U7
 */
static void mm_CoppersmithWinograd(float* matA, float* matB, float* matC, const int M, const int N, const int K,
				   const int strideA, const int strideB, const int strideC) {
  // LEAF_SIZE = 2 OR DIMENSIONS OF MATRIX ARE ODD
  if ((M <= 2) || (M % 2 != 0 || N % 2 != 0 || K % 2 != 0)) {
    return mm_basic(matA, matB, matC, M, N, K, strideA, strideB, strideC);
  }

  float* S1 = (float*)malloc((M / 2) * (K / 2) * sizeof(float));
  float* S2 = (float*)malloc((M / 2) * (K / 2) * sizeof(float));
  float* S3 = (float*)malloc((M / 2) * (K / 2) * sizeof(float));
  float* S4 = (float*)malloc((M / 2) * (K / 2) * sizeof(float));
  {
    for (int i = 0; i < M / 2; i++) {
      for (int j = 0; j < K / 2; j++) {
	int idxA, offset, idxS = i * (K / 2) + j;

	//S1     = A21 + A22
	idxA = (i + (M / 2)) * strideA + j;
	offset = K / 2;
	S1[idxS] = matA[idxA] + matA[idxA + offset];

	//S2     = S1 - A11
	idxA = i * strideA + j;
	S2[idxS] = S1[idxS] - matA[idxA];

	//S3     = A11 - A21
	offset = (M / 2) * strideA;
	S3[idxS] = matA[idxA] - matA[idxA + offset];

	//S4     = A12 - S2
	idxA = i * strideA + (K / 2) + j;
	S4[idxS] = matA[idxA] - S2[idxS];
      }
    }
  }

  float* T1 = (float*)malloc((K / 2) * (N / 2) * sizeof(float));
  float* T2 = (float*)malloc((K / 2) * (N / 2) * sizeof(float));
  float* T3 = (float*)malloc((K / 2) * (N / 2) * sizeof(float));
  float* T4 = (float*)malloc((K / 2) * (N / 2) * sizeof(float));
  {
    for (int i = 0; i < K / 2; i++) {
      for (int j = 0; j < N / 2; j++) {
	int idxB, offset, idxT = i * (N / 2) + j;

	//T1     = B12 - B11
	idxB = i * strideB + j;
	offset = (N / 2);
	T1[idxT] = matB[idxB + offset] - matB[idxB];

	//T2     = B22 - T1
	idxB = (i + (K / 2)) * strideB + (N / 2) + j;
	T2[idxT] = matB[idxB] - T1[idxT];

	//T3     = B22 - B12
	idxB = i * strideB + (N / 2) + j;
	offset = ((K / 2)) * strideB;
	T3[idxT] = matB[idxB + offset] - matB[idxB];

	//T4     = T2 - B21
	idxB = (i + (K / 2)) * strideB + j;
	T4[idxT] = T2[idxT] - matB[idxB];
      }
    }
  }

  //M1 = A11 * B11
  float* M1 = (float*)malloc((M / 2) * (N / 2) * sizeof(float));
  mm_CoppersmithWinograd(matA, matB, &M1[0], M / 2, N / 2, K / 2, strideA, strideB, N / 2);

  //M2 = A12 * B21
  float* M2 = (float*)malloc((M / 2) * (N / 2) * sizeof(float));
  mm_CoppersmithWinograd(&matA[K / 2], &matB[(K / 2) * strideB], &M2[0], M / 2, N / 2, K / 2, strideA, strideB, N / 2);

  //M3 = S4 * B22
  float* M3 = (float*)malloc((M / 2) * (N / 2) * sizeof(float));
  mm_CoppersmithWinograd(&S4[0], &matB[(K / 2) * strideB + (N / 2)], &M3[0], M / 2, N / 2, K / 2, K / 2, strideB, N / 2);

  //M4 = A22 * T4
  float* M4 = (float*)malloc((M / 2) * (N / 2) * sizeof(float));
  mm_CoppersmithWinograd(&matA[(M / 2) * strideA + (K / 2)], &T4[0], &M4[0], M / 2, N / 2, K / 2, strideA, N / 2, N / 2);

  //M5 = S1 * T1
  float* M5 = (float*)malloc((M / 2) * (N / 2) * sizeof(float));
  mm_CoppersmithWinograd(&S1[0], &T1[0], &M5[0], M / 2, N / 2, K / 2, K / 2, N / 2, N / 2);

  //M6 = S2 * T2
  float* M6 = (float*)malloc((M / 2) * (N / 2) * sizeof(float));
  mm_CoppersmithWinograd(&S2[0], &T2[0], &M6[0], M / 2, N / 2, K / 2, K / 2, N / 2, N / 2);

  //M7 = S3 * T3
  float* M7 = (float*)malloc((M / 2) * (N / 2) * sizeof(float));
  mm_CoppersmithWinograd(&S3[0], &T3[0], &M7[0], M / 2, N / 2, K / 2, K / 2, N / 2, N / 2);

  //C11 = U1 = M1 + M2
  //C12 = U5 = U4 + M3 = U2 + M5 + M3 = M1 + M6 + M5 + M3
  //C21 = U6 = U3 - M4 = U2 + M7 - M4 = M1 + M6 + M7 - M4
  //C22 = U7 = U3 + M5 = U2 + M7 + M5 = M1 + M6 + M7 + M5
  for (int i = 0; i < M / 2; i++) {
    for (int j = 0; j < N / 2; j++) {
      int idx = i * (N / 2) + j;
      matC[i * strideC + j] = M1[idx] + M2[idx];
      matC[i * strideC + j + (N / 2)] = M1[idx] + M6[idx] + M5[idx] + M3[idx];
      matC[(i + (M / 2)) * strideC + j] = M1[idx] + M6[idx] + M7[idx] - M4[idx];
      matC[(i + (M / 2)) * strideC + j + (N / 2)] = M1[idx] + M6[idx] + M7[idx] + M5[idx];
    }
  }
  free(S1);           S1 = NULL;
  free(S2);           S2 = NULL;
  free(S3);           S3 = NULL;
  free(S4);           S4 = NULL;
  free(T1);           T1 = NULL;
  free(T2);           T2 = NULL;
  free(T3);           T3 = NULL;
  free(T4);           T4 = NULL;
  free(M1);           M1 = NULL;
  free(M2);           M2 = NULL;
  free(M3);           M3 = NULL;
  free(M4);           M4 = NULL;
  free(M5);           M5 = NULL;
  free(M6);           M6 = NULL;
  free(M7);           M7 = NULL;
}

static void mm_CoppersmithWinograd(float* matA, float* matB, float* matC, const int M, const int N, const int K) {
  mm_CoppersmithWinograd(matA, matB, matC, M, N, K, K, N, N);
}


static void mm_test(int M, int N, int K, int rangeTop) {
  unsigned seed = time(0);
  srand(seed);
  clock_t start, end;
  for (int i = 0; i < 1; i++) {
    float* mA = (float*)malloc(M * K * sizeof(float));
    float* mB = (float*)malloc(K * N * sizeof(float));
    float* mE = (float*)malloc(M * N * sizeof(float));
    for (int j = 0; j < M * K; j++) {
      mA[j] = rand() % rangeTop;
    }
    for (int j = 0; j < K * N; j++) {
      mB[j] = rand() % rangeTop;
    }
    
    start = clock();
    mm_CoppersmithWinograd(mA, mB, mE, M, N, K);
    end = clock();
    auto endtime = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Winograd%d time: %fms\n", i, endtime * 1000);
    printf("\n");
    free(mA);
    free(mB);
    free(mE);
  }
}

int main() {
  int M, N, K, rangeTop;
  M = 200;
  N = 200;
  K = 200;
  rangeTop = 100;
  printf("Start\n");
  mm_test(M, N, K, rangeTop);
  return 0;
}
