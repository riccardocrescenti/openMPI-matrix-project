#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>

#define DIMBUFF 100000

//Allocation of a memory contiguous 2D array of float numbers
float **alloc(int rows, int cols) {
    float *data = (float *)malloc(rows*cols*sizeof(float));
    float **array= (float **)malloc(rows*sizeof(float*));
    int i;
    for (i=0; i<rows; i++) {
    	array[i] = &(data[cols*i]);
    }

    return array;
}

//Deallocation of a memory contiguous 2D array of float numbers
void free2D(float **mat) {
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
	matrix = alloc(rows,columns);
	
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

//Function which returns the product of two given matrices
float **matrixProduct(float **matrix1, int nrows1, int ncolumns1, float **matrix2, int nrows2, int ncolumns2) {
	//Check matrices dimensions to find if a multiplication in possible
	if(ncolumns1 != nrows2) {
		puts("It is not possible to multiply these matrices");
		return NULL;
	}
	
	float **ris;
	
	//Matrix allocation, first rows the columns
	ris = alloc(nrows1, ncolumns2);
	
	//Computing the product
	for(int i = 0; i < nrows1; i++) {
		for(int j = 0; j < ncolumns2; j++) {
			float sum = 0;
			for(int k = 0; k < ncolumns1; k++) {
				sum += matrix1[i][k] * matrix2[k][j];
			}
			ris[i][j] = sum;
		}
	}
	
	return ris;
}




void main(int argc, char *argv[]) {
	//Checking the arguments
	if(argc != 3) {
		puts("You have to pass two files which contain the matrices");
		return;
	}
	
	//Declaration of MPI variables
	MPI_Status status;
	int myRank, P;
	int q, r;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);

	//Declaration of variables
	float **matrix1, **matrix2;
	float **recvMatrix1;
	float **localResult;
	float **resultMatrix;
	int r1, c1, r2, c2;
	
	//process 0 reads matrices from files
	if (myRank == 0){
		matrix1 = readMatrixFromFile(argv[1], &r1, &c1);
		if(matrix1 == NULL) {
			printf("File %s is not correct\n", argv[1]);
			return;
		}
		
		matrix2 = readMatrixFromFile(argv[2], &r2, &c2);
		if(matrix2 == NULL) {
			printf("File %s is not correct\n", argv[2]);
			return;
		}	
		
		//computing size of splitted matrix
		q = r1 / P;
		r = r1 % P;
		
		if (q < 1) {
			q = 1;
			r = 0;
		}
		
		resultMatrix = alloc(r1, c2);
		
		//Printing the two matrices
		printMatrix(matrix1, r1, c1);
		printMatrix(matrix2, r2, c2);
	}
	
	
	
	//broadcasting dimensions
	MPI_Bcast(&r1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&c1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&r2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&c2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&q, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&r, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	//Matrices allocation
	if (myRank != 0){
		matrix1 = alloc(r1, c1);
		matrix2 = alloc(r2, c2);
	}
	recvMatrix1 = alloc(q, c1);
	localResult = alloc(q, c2);
	
	//Sending matrices to all processes
	MPI_Bcast(&matrix2[0][0], r2 * c2, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Scatter(matrix1[0], q * c1, MPI_FLOAT, &recvMatrix1[0][0], q * c1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
	//Computing local products
	localResult = matrixProduct(recvMatrix1, q, c1, matrix2, r2, c2);

	//Sending local products back to master process
	MPI_Gather(localResult[0], q * c2, MPI_FLOAT, &resultMatrix[0][0], q * c2, MPI_FLOAT, 0, MPI_COMM_WORLD); 

	//Master process compute the product for last lines (in case matrix dimension was not a multiple of the number of processes) and prints the result
	if (myRank == 0){
		if (r != 0){
			float **lastLines;
			float **lastResults;
			
			lastLines = alloc(r, c1);
			lastResults = alloc(r, c2);
			
			int startRows = r1 - r;
			
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < c1; j++){
					lastLines[i][j] = matrix1[startRows][j];
				}
				startRows++;
			}
			
			lastResults = matrixProduct(lastLines, r, c1, matrix2, r2, c2);
	
			startRows = r1 - r;
			
			for (int i = 0; i < r ; i++){
				for (int j = 0; j < c2 ; j++){
					resultMatrix[startRows][j] = lastResults[i][j];
				}
				startRows++;
			}			
			
			printMatrix(resultMatrix, r1, c2);
			free2D(lastLines);
			free2D(lastResults);
			free2D(matrix1);
			free2D(resultMatrix);
		} else {
		
			printMatrix(resultMatrix, r1, c2);
			free2D(matrix1);
			free2D(resultMatrix);
		}		
	}

	free2D(matrix2);
	free2D(recvMatrix1);
	free2D(localResult);

	MPI_Finalize();
	
}
