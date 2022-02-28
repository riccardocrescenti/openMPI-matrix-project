#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define NAME 30

//Creation of a matrix with size rows x columns, stored in a file .txt
int main() {
	//Declaration of variables
	int rows, columns;
	char fileName[NAME];
	
	//Seed for casual numbers placed in the matrix
	srand(time(NULL));
	
	//Setting matrix dimensions and file name
	puts("Insert a number of rows:");
	scanf("%d", &rows);
	while(rows <= 0) {
		puts("The number of rows must be positive, insert here:");
		scanf("%d", &rows);
	}
	
	puts("Insert a number of columns:");
	scanf("%d", &columns);
	while(columns <= 0) {
		puts("The number of columns must be positive, insert here:");
		scanf("%d", &columns);
	}
	
	puts("Insert a name for the file:");
	scanf("%s", fileName);
	while(fileName == NULL || strlen(fileName) >= 30) {
		puts("This name is not valid, insert a new one here:");
		scanf("%s", fileName);
	}
	
	//Writing into the file:
	//First line contains the number of rows and columns
	//The rest of the file is a columns with a single value in each row
	FILE *f = fopen(fileName, "w");
	
	fprintf(f, "%d %d\n", rows, columns);
	for(int i = 0; i < rows * columns; i++) {
		float casualNumber = (rand() % 10) + 1;
		fprintf(f, "%f\n", casualNumber);
	}
	
	fclose(f);

	return 0;
}
