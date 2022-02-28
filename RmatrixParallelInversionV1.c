#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#define DIMBUFF 1000000
#define TOL 0.0001

//Allocation of a memory contiguous array of int
int *allocVector(int dim){
	int *array = (int *)malloc(dim * sizeof(int));
	
	return array;
}

//Allocation of a memory contiguous array of float
float *allocVectorFloat(int dim){
	float *array = (float *)malloc(dim * sizeof(float));
	
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

//Deallocation of a memory contiguous array of float
void free1DFloat(float *vector){
	free(vector);
}

//Deallocation of a memory contiguous 2D array of float
void free2D(float **mat){
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
	
	//Matrix allocation
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
	
	//Declaration of MPI variables
	MPI_Status status;
	int myRank, size, retVal;
	int q, r;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int n, columns;
	float **A, **Ainv;
	int *P;
	float *colArray, *diagArray;
	float **matrixPortion;
	int tol = TOL;
	
	//Master process reads the matrix from file
	if (myRank == 0){
		A = readMatrixFromFile(argv[1], &n, &columns);
		if(A == NULL) {
			printf("File %s is not correct\n", argv[1]);
			MPI_Finalize();
			return 0;
		}

		//Checking if A is a square matrix
		if(n != columns) {
			puts("This is not a square matrix");
			MPI_Finalize();
			return 0;
		}
		
		P = allocVector(n + 1);
		Ainv = alloc(n, n);
		
		//Initialization of P vector
		for (int i = 0; i <= n; i++){
			P[i] = i;
		}
		
		printMatrix(A, n, n);	
	}

	//Broadcasting matrix dimension
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Computation of LU decomposition (there are loop carried dependencies so we cannot parallelize the entire process)
	for (int col = 0; col < n - 1; col++) {
		//Matrix scattering among processes for pivot computation
		if (myRank == 0) {
			q = (n - col) / size;
			r = (n - col) % size;

			colArray = allocVectorFloat(q + r);
			
			int newIndex = col;
			
			for (int i = 0; i < q + r; i++){
				colArray[i] = A[newIndex][col];
				newIndex++;
			}
			
			for (int p = 1; p < size; p++){
				for (int i = q * p + r + col; i < q * p + r + col + q; i++){
					retVal = MPI_Send(&A[i][col], 1, MPI_FLOAT, p, 555, MPI_COMM_WORLD);
				}
			}
		} else {
			q = (n - col) / size;
			r = (n - col) % size;	
			colArray = allocVectorFloat(q);	
		
			for (int i = 0; i < q; i++){		
				retVal = MPI_Recv(&colArray[i], 1, MPI_FLOAT, 0, 555, MPI_COMM_WORLD, &status);
			}
		}
		
		//Computation of column maximum and gathering to master process
		int startIndex;
		float maxA, absA;
		int imax, j;
		float *ptr;
		int indexVector[size];
		
		for (int i = 0; i < size; i++){
			indexVector[i] = -1;		
		}
		
		if (myRank == 0){
			startIndex = col;
			maxA = 0;
			imax = startIndex;
			for (int i = 0; i < q + r; i++){
				if ((absA = fabs(colArray[i])) > maxA){
					maxA = absA;
					imax = i;
				}
			}
			imax += startIndex;	
		} else {
			startIndex = q * myRank + r + col;
			maxA = 0;
			imax = startIndex;
			for (int i = 0; i < q; i++){
				if ((absA = fabs(colArray[i])) > maxA){
					maxA = absA;
					imax = i;
				}
			}
			imax += startIndex;	
		}
		//Gathering only if other processes have done computations
		if (col <= n - size){
			MPI_Gather(&imax, 1, MPI_INT, &indexVector[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (myRank == 0){
				for (int i = 0; i < size; i++){
					if (A[indexVector[i]][col] > maxA){
						maxA = A[indexVector[i]][col];
						imax = indexVector[i];
					}
				}
			}
		}
		
		//Master process switch rows
		if (myRank == 0){
			if (maxA < tol){
				puts("Matrice degenerata");
				MPI_Finalize();
				return 0;
			}
			
			if (imax != col) {
		        //pivoting P
				j = P[col];
				P[col] = P[imax];
				P[imax] = j;

				//pivoting rows of A
		        ptr = A[col];
		        A[col] = A[imax];
		        A[imax] = ptr;
		        //counting pivots starting from N (for determinant)
		        P[n]++;
			}
		}
			
		
		free1DFloat(colArray);	
		
		//Starting LU DECOMPOSITION
		float mainRow[n];
		
		//Broadcasting the main row and scattering the rest of the matrix
		if (myRank == 0){
			for (int i = 0; i < n; i++){
				mainRow[i] = A[col][i];
			}
		}
		MPI_Bcast(&mainRow, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		if (myRank == 0) {
			// -1 to skip the main row already broadcasted
			q = (n - col - 1) / size;
			r = (n - col - 1) % size;
			
			matrixPortion = alloc(q + r, n);
			
			int nuovoIndex = col + 1;
			for (int i = 0; i < q + r; i++){
				matrixPortion[i] = A[nuovoIndex];
				nuovoIndex++;
			}
			
			
			for (int p = 1; p < size; p++){
				for (int i = q * p + r + col + 1; i < q * p + r + col + 1 + q; i++){
					retVal = MPI_Send(&A[i][0], n, MPI_FLOAT, p, 555, MPI_COMM_WORLD);
				}
			}
		} else {
			q = (n - col - 1) / size;
			r = (n - col - 1) % size;	
			matrixPortion = alloc(q, n);	
		
			for (int i = 0; i < q; i++){		
				retVal = MPI_Recv(&matrixPortion[i][0], n, MPI_FLOAT, 0, 555, MPI_COMM_WORLD, &status);
			}
			
		}
		
		//Computing l coefficient and updating rows
		if (myRank == 0){
			for (int i = 0; i < q + r; i++){
				matrixPortion[i][col] /= mainRow[col];
				
				for (int k = col + 1; k < n; k++) {
					matrixPortion[i][k] -= matrixPortion[i][col] * mainRow[k];
				}
				
			}
		} else {
			for (int i = 0; i < q; i++){
				matrixPortion[i][col] /= mainRow[col];
				
				for (int k = col + 1; k < n; k++) {
					matrixPortion[i][k] -= matrixPortion[i][col] * mainRow[k];
				}
			}
		}
		
		startIndex = col + 1;
		
		//Send matrix portions back to the master
		if (myRank == 0){
			for (int i = 0; i < q + r; i++){
				for (int j = 0; j < n; j++){
					A[startIndex][j] = matrixPortion[i][j];
				}
				startIndex++;
			}
		}

		if (col < n - size){
			if (myRank != 0){
				for (int i = 0; i < q; i++){
					retVal = MPI_Send(&matrixPortion[i][0], n, MPI_FLOAT, 0, 555, MPI_COMM_WORLD);
				}
			}
			if (myRank == 0){
				for (int p = 1; p < size; p++) {
					for (int i = 0; i < q; i++){
						retVal = MPI_Recv(&A[col + 1 + r + p * q + i][0], n, MPI_FLOAT, p, 555, MPI_COMM_WORLD, &status);
					}
				}
			}
		}
	}

	//End of LU Decomposition
	
	
	//Check if matrix A is singular computing the determinant (if true the program stops because A is not invertible)	
	if (myRank == 0){
		q = n / size;
		r = n % size;
		
		diagArray = allocVectorFloat(q + r);
		
		for (int i = 0; i < q + r; i++){
			diagArray[i] = A[i][i];
		}		
		for (int p = 1; p < size; p++){
			for (int i = q * p + r; i < q * p + r + q; i++){
				retVal = MPI_Send(&A[i][i], 1, MPI_FLOAT, p, 555, MPI_COMM_WORLD);
			}
		}		
	} else {
		q = n / size;
		r = n % size;
		diagArray = allocVectorFloat(q);
	
		for (int i = 0; i < q; i++){		
			retVal = MPI_Recv(&diagArray[i], 1, MPI_FLOAT, 0, 555, MPI_COMM_WORLD, &status);
		}
	}
	
	float det = diagArray[0];
	float diagProduct[size];
	
	//Computing multiplication of diagonal elements and sending local results back to the master process
	if (myRank == 0){
		for (int i = 1; i < q + r; i++){
			det *= diagArray[i];
		}
	}
	else {
		for (int i = 1; i < q; i++){
			det *= diagArray[i];			
		}
	}
	if (n >= size){
		MPI_Gather(&det, 1, MPI_FLOAT, &diagProduct[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (myRank == 0){
			for (int i = 1; i < size; i++){
				det *= diagProduct[i];
			}
			det = (P[n] - n) % 2 == 0 ? det : -det;
		}
	} else {
		if (myRank == 0){
			det = (P[n] - n) % 2 == 0 ? det : -det;
		}
	}
	free1DFloat(diagArray);	
	
	if (myRank == 0){
		if (det == 0){
			puts("La matrice non è invertibile");
			MPI_Finalize();
			return 0;
		}

	}
	//End of determinant checking

	//Inversion
	//There are dependencies between rows but not between columns, so we split matrix by columns
	
	//Broadcasting of A' and P
	if (myRank != 0){
		A = alloc(n, n);
		P = allocVector(n + 1);
	}
		
	if (myRank == 0){
		for (int p = 1; p < size; p++){
			for (int i = 0; i < n; i++){
				for (int j = 0; j < n; j++){
					retVal = MPI_Send(&A[i][j], 1, MPI_FLOAT, p, 555, MPI_COMM_WORLD);
				}
			}
		}
	}
	
	if (myRank != 0){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				retVal = MPI_Recv(&A[i][j], 1, MPI_FLOAT, 0, 555, MPI_COMM_WORLD, &status);
			}
		}

	}

	MPI_Bcast(&P[0], n + 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	float **recvMat;
	
	q = n / size;
	r = n % size;
	
	if (myRank == 0){
		recvMat = alloc(n, q + r); 	
	} else {
		recvMat = alloc(n, q);
	}
	//Initializing and sending matrix divided by columns, then computing inversion
	if (myRank == 0){
		for (int j = 0; j < q + r; j++){
			for (int i = 0; i < n; i++){
				recvMat[i][j] = P[i] == j ? 1.0 : 0.0;

				for (int k = 0; k < i; k++){
					recvMat[i][j] -= A[i][k] * recvMat[k][j];
				}
			}
			
			for (int i = n - 1; i >= 0; i--){
				for (int k = i + 1; k < n; k++){
					recvMat[i][j] -= A[i][k] * recvMat[k][j];
				}
				recvMat[i][j] /= A[i][i];
			}
		}
	}
	if (n >= size){
	
		int realColumn = q * myRank + r;
	
		if (myRank != 0){
			for (int j = 0; j < q; j++){
				for (int i = 0; i < n; i++){
					recvMat[i][j] = P[i] == realColumn ? 1.0 : 0.0;

					for (int k = 0; k < i; k++){
						recvMat[i][j] -= A[i][k] * recvMat[k][j];
					}
				}
				
				for (int i = n - 1; i >= 0; i--){
					for (int k = i + 1; k < n; k++){
						recvMat[i][j] -= A[i][k] * recvMat[k][j];
					}
					recvMat[i][j] /= A[i][i];
				}
				realColumn++;
			}
		}
	}
	
	//Gathering local results back to the master
	if (myRank == 0){
		for (int j = 0; j < q + r; j++){
			for (int i = 0; i < n; i++){
				Ainv[i][j] = recvMat[i][j];
			}
		}
	} 
	if (myRank != 0) {
		for (int j = 0; j < q; j++){
			for (int i = 0; i < n; i++){
				retVal = MPI_Send(&recvMat[i][j], 1, MPI_FLOAT, 0, 555, MPI_COMM_WORLD);
			}
		}
	}
	
	if (myRank == 0){
		for (int p = 1; p < size; p++){
			for (int j = q * p + r; j < q * p + r + q; j++){
				for (int i = 0; i < n; i++){
					retVal = MPI_Recv(&Ainv[i][j], 1, MPI_FLOAT, p, 555, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	

	//Master process prints the result
	if (myRank == 0){
		printMatrix(Ainv, n, n);
		free1D(P);
		free2D(Ainv);
	}
	
	if (myRank != 0){
		free2D(A);
	}
	MPI_Finalize();
	
	return 0;
}
