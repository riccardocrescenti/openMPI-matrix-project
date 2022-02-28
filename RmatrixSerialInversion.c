#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define DIMBUFF 100000

//Computation of the determinant of matrix A by multiplying all the diagonal elements of the LU decomposed matrix
float LUPDeterminant(float **A, int *P, int N) {

    float det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    return (P[N] - N) % 2 == 0 ? det : -det;
}

//Algorithm for matrix inversion starting from the LU decomposed matrix
void LUPInvert(float **A, int *P, int N, float **IA) {
  
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

//LU decomposition of matrix A (executed in place)
int LUPDecompose(float **A, int N, float Tol, int *P) {

    int i, j, k, imax; 
    float maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) { 
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
        
    }

    return 1;  //decomposition done 
}

//Allocation of a memory contiguous array of int
int *allocVector(int dim){
	int *array = (int *)malloc(dim * sizeof(int));
	
	return array;
}

//Allocation of a memory contiguous 2D array of float
float **alloc(int rows, int cols)
{
    float *data = (float *)malloc(rows*cols*sizeof(float));
    float **array= (float **)malloc(rows*sizeof(float*));
    int i;
    for (i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

//Deallocation of a memory contiguous array of int
void free1D(int *vector){
	free(vector);
}

//Deallocation of a memory contiguous 2D array of float
void free2D(float **mat)
{
    free(mat[0]);
    free(mat);
}

//Function which reads a matrix from a file and return it
float **readMatrixFromFile(char *path, int *nrows, int *ncolumns) {
	float **matrix;
	FILE *f;
	char buffer[DIMBUFF];
	int rows, columns;
	
	//File opening
	if(!(f = fopen(path, "r"))) {
		puts("File not found");
		return NULL;
	}
	
	//Reading the first line and passing arguments to the main program
	fgets(buffer, sizeof(buffer), f);
	sscanf(buffer, "%d %d", &rows, &columns);
	*nrows = rows;
	*ncolumns = columns;
	
	//Matrix allocation, first rows the columns
	matrix = alloc(rows, columns);
	
	//Reading matrix values
	int i = 0, j = 0;
	int counter = 0;
	while(fgets(buffer, sizeof(buffer), f)) {
		sscanf(buffer, "%f", &matrix[i][j]);
		j++;
		counter++;
		
		//When j == columns a row is filled, so I reset j and read the following row
		if(j == columns) {
			j = 0;
			i++;
		}
	}
	
	fclose(f);
	
	//Check validity and return the result
	if(counter != rows * columns) {
		return NULL;
	}
	return matrix;
}

//Function which prints a matrix and put a \n at the end
void printMatrix(float **matrix, int nrows, int ncolumns) {
	for(int i = 0; i < nrows; i++) {
		for(int j = 0; j < ncolumns; j++) {
			printf("%.3f ", matrix[i][j]);
		}
		puts("");
	}
	puts("");
}

//Inversion of a matrix NxN using LU decomposition method
int main(int argc, char *argv[]) {
	//Checking the arguments
	if(argc != 2) {
		puts("You have to pass a file which contains the matrix");
		return -1;
	}
	
	int n, columns;
	float **A, **Ainv;
	int *P;

	//Reading matrix from file
	A = readMatrixFromFile(argv[1], &n, &columns);
	if(A == NULL) {
		printf("File %s is not correct\n", argv[1]);
		return -1;
	}

	//Checking if A is a square matrix
	if(n != columns) {
		puts("This is not a square matrix");
		return -1;
	}
	
	P = allocVector(n + 1);
	Ainv = alloc(n, n);
	printMatrix(A, n, n);
	
	//Computing LU decomposition of matrix A
	LUPDecompose(&A[0], n, 0.0001, &P[0]);
	
//	puts("LU Decomposition");
//	printMatrix(A, n, n);

	float determinant;
	
	determinant = LUPDeterminant(A, P, n);
	
	//Checking if matrix A is singular (in that case the program stops because A is not invertible)
	if (determinant == 0){
		puts("La matrice non è invertibile");
		return -1;
	}

	//Computing matrix inversion
	LUPInvert(A, P, n, &Ainv[0]);
	
	//Printing the result
	printMatrix(Ainv, n, n);
	
	
	free1D(P);
	free2D(Ainv);

	return 0;
}
