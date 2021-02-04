/*  Name: Nick Dreitzler
    Username: ndreitz
    CUID: C15564443
    Class: ECE 4730
    Professor: Walt Ligon
    Semester: Fall2020
    Filename: floyd-serial.c

    Purpose: Reads 2 matrix files, multiplies them and writes the
        output to a binary file 

    usage: ./mm-serial <input_file1> <input_file2> <output_file>
    
    Output: result of matrix multiplication
 */

#include "matrix.h"

int check_dims(char *file_name1, char *file_name2);
void matrix_mult( my_matrix *Matrix1, my_matrix *Matrix2, my_matrix *Result);

int main(int argc, char **argv)
{
    my_matrix *Matrix1 = (my_matrix *)malloc(sizeof(my_matrix));
    my_matrix *Matrix2 = (my_matrix *)malloc(sizeof(my_matrix));
    my_matrix *Result = (my_matrix *)malloc(sizeof(my_matrix));
    clock_t exe_time;

    if(argc != 4)
    {
        printf("Usage:\n\t./floyd-serial <input file1> <input_file2> <output file>\n");
        exit(1);
    }

    if(!check_dims(argv[1], argv[2]))
    {
        printf("Number of columns in matrix 1 and number of rows in matrix 2 must match\n");
        exit(0);
    }

    read_matrix(argv[1], Matrix1);
    read_matrix(argv[2], Matrix2);

    exe_time = -clock();
    matrix_mult(Matrix1, Matrix2, Result);
    exe_time += clock();

    write_matrix(argv[3], Result);

    printf("%5s\t%5s\t%5s\t%13s\n", "","size", "n", "exe_time (ms)");
    printf("%5s\t%5d\t%5d\t%13f\n", "data",1 , Result->dims[ROW], ((double)exe_time*1000)/CLOCKS_PER_SEC);

    free_matrix(&Matrix1);
    free_matrix(&Matrix2);
    free_matrix(&Result);
}

/*Checks that the dimensions of the 2 matrices are a match for multiplication.
    Number of columns in matrix 1 must equal the number of rows in matrix 2.
    returns 1 if they match 2 if they don't.*/
//int check_dims(char *file_name1, char *file_name2, FILE **M1_fptr, FILE **M2_fptr)
int check_dims(char *file_name1, char *file_name2)
{
    FILE *M1_fptr;
    FILE *M2_fptr;
    int M1_dims[2];
    int M2_dims[2];

    OPEN_READ(M1_fptr, file_name1);
    OPEN_READ(M2_fptr, file_name2);

    fread(M1_dims, 2, sizeof(int), M1_fptr);
    fread(M2_dims, 2, sizeof(int), M2_fptr);

    // printf("M1 %d x %d\n", M1_dims[ROW], M1_dims[COL]);
    // printf("M2 %d x %d\n", M2_dims[ROW], M2_dims[COL]);
    // printf("%d\n", (M1_dims[COL] == M2_dims[ROW]));

    CLOSE(M1_fptr);
    CLOSE(M2_fptr);

    return (M1_dims[COL] == M2_dims[ROW]);
}

/*Multiplies two matrices together and stores the result in a new matrix
*/
void matrix_mult(
    my_matrix *Matrix1, /*IN  - multiplicand*/
    my_matrix *Matrix2, /*IN  - multiplier*/
    my_matrix *Result)  /*OUT - result of multiplication*/
{
    int i,j,k;

    Result->dims[ROW] = Matrix1->dims[ROW];
    Result->dims[COL] = Matrix2->dims[COL]; 

    Result->block_dims[ROW] = Matrix1->block_dims[ROW];
    Result->block_dims[COL] = Matrix2->block_dims[COL]; 

    calloc_matrix(Result);

    for (i = 0; i < Matrix1->dims[ROW]; ++i)
        for (j = 0; j < Matrix2->dims[COL]; ++j)
            for (k = 0; k < Matrix1->dims[COL]; ++k) 
                Result->data[i][j] += Matrix1->data[i][k] * Matrix2->data[k][j];


//     *result_storage = storage;
//     *result = result;
}
