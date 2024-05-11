*********************************************
* Joey Mulé, Chris Bispels
* University of Maryland, Baltimore County
* CMSC 441-Algorithms, Dr. Sherman	  
* 5/10/24		   
*********************************************
-----------------------------------------------
- BAM vs SAM vs CWA                           
-----------------------------------------------
- Algorithm Running Time and Space comparison 
-----------------------------------------------

The following represents our code for this project:
------------------------------------------------------------------------------------------------------------------------------------
- BAM_SAM_CWA.cpp: Combined file which includes all three algorithms running 5 total trials (on each), calculating		   
- running time via timestamping. Each matrix is randomized per trial. The order of which the algorithms are ran are: BAM, SAM, CWA 
------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------
- BAMspace.cpp: File containing ONLY the BAM algorithm, which is set to run only once, and is intended to be ran with valgrind to 
- capture space usage.														  
- After Compilation:												                  
- valgrind ./BAMspace												                  
-----------------------------------------------------------------------------------------------------------------------------------
- SAMspace.cpp: File containing ONLY the SAM algorithm, which is set to run only once, and is intended to be ran with valgrind to 
- capture space usage.											 	                  
- After Compilation:													          
- valgrind ./SAMspace														  
-----------------------------------------------------------------------------------------------------------------------------------
- CWAspace.cpp: File containing ONLY the CWA algorithm, which is set to run only once, and is intended to be ran with valgrind to 
- capture space usage.														  
- After Compilation:													 	  
- valgrind ./CWAspace														  
-----------------------------------------------------------------------------------------------------------------------------------

  TO CHANGE THE MATRIX SIZE: 
==============================

Inside of the main() func

```
  int main(){;
    int M, N, K, rangeTop;
    // MxK * KxN
    M = 200; <------------------- CHANGE
    N = 200; <------------------------- THESE
    K = 200; <--------------------------------- VALUES
	// matX[i][j] values - (0 - 100)
    rangeTop = 100;
	printf("start algorithms\n\n");
    mm_test(M, N, K, rangeTop);
    return 0;
}
```


