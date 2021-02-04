/*  Name: Nick Dreitzler
    Username: ndreitz
    CUID: C15564443
    Class: ECE 4730
    Professor: Walt Ligon
    Semester: Fall2020
    Filename: print-matrix.c

    Purpose: To read a binary file of a r by c matrix .
        Data are doubles and the ints r anc c are stored as 
        the first and second words of the file.

    Usage: ./print-matrix <input file>
    
    Output: Print out of given matrix
 */

#include "matrix.h"

int main(int argc, char **argv)
{
    my_matrix *Matrix = (my_matrix *)malloc(sizeof(my_matrix));
    // dtype *storage;
    // dtype **matrix;
    // int Matrix->dims[2];

    if(argc != 2)
    {
        printf("Usage:\n\t./print-file <input file>\n");
        exit(1);
    }

    //Read adjacency matrix from argv[1]
    read_matrix(argv[1], Matrix);

    print_matrix(Matrix);

    //Clean up
    // FREE(storage);
    // FREE(matrix);
    free_matrix(&Matrix);
    return(0);
}