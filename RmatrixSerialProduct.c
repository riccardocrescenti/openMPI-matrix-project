#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define DIMBUFF 100000

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
	matrix = (float **) malloc(rows * sizeof(float *));
	if(matrix == NULL) {
		puts("Error in matrix allocation");
		return NULL;
	}
	
	for(int i = 0; i < rows; i++) {
		matrix[i] = (float *) malloc(columns * sizeof(float));
		if(matrix[i] == NULL) {
			printf("Error in allocation of row %d\n", i);
			return NULL;
		}
	}
	
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

//Function which frees the memory allocated before for the given matrix
void freeMemory(float **matrix, int nrows) {
	for(int i = 0; i < nrows; i++) {
		free(matrix[i]);
	}
	free(matrix);
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
	ris = (float **) malloc(nrows1 * sizeof(float *));
	if(ris == NULL) {
		puts("Error in matrix allocation");
		return NULL;
	}
	
	for(int i = 0; i < nrows1; i++) {
		ris[i] = (float *) malloc(ncolumns2 * sizeof(float));
		if(ris[i] == NULL) {
			printf("Error in allocation of row %d\n", i);
			return NULL;
		}
	}
	
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

//Product of two matrices with dimensions NxM and MxQ --> Result dimensions are NxQ
int main(int argc, char *argv[]) {
	//Checking the arguments
	if(argc != 3) {
		puts("You have to pass two files which contain the matrices");
		return -1;
	}
	
	//Declaration of variables
	float **matrix1, **matrix2;
	float **resultMatrix;
	int rows1, columns1, rows2, columns2;
	
	//Reading matrices from files
	matrix1 = readMatrixFromFile(argv[1], &rows1, &columns1);
	if(matrix1 == NULL) {
		printf("File %s is not correct\n", argv[1]);
		return -1;
	}
	
	matrix2 = readMatrixFromFile(argv[2], &rows2, &columns2);
	if(matrix2 == NULL) {
		printf("File %s is not correct\n", argv[2]);
		return -1;
	}
	
	//Printing the two matrices
	printMatrix(matrix1, rows1, columns1);
	printMatrix(matrix2, rows2, columns2);
	
	//Computing the product
	resultMatrix = matrixProduct(matrix1, rows1, columns1, matrix2, rows2, columns2);
	
	if(resultMatrix == NULL) {
		return -1;
	}
	
	//Printing the result
	printMatrix(resultMatrix, rows1, columns2);
	
	//Freeing memory allocated in file read and product computation
	freeMemory(matrix1, rows1);
	freeMemory(matrix2, rows2);
	freeMemory(resultMatrix, rows1);

	return 0;
}
