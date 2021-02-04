/*  Name: Nick Dreitzler
    Username: ndreitz
    CUID: C15564443
    Class: ECE 4730
    Professor: Walt Ligon
    Semester: Fall2020
    Filename: graph.h
 */

//Guard
#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#define ROW 0
#define COL 1
#define LOWER 0
#define UPPER 1
#define NDIMS 2

#define FREE(x)                                         \
do {                                                    \
    free(x);                                            \
    x = NULL;                                           \
} while(0)

#define CLOSE(x)                                        \
do {                                                    \
    fclose(x);                                          \
    x = NULL;                                           \
} while(0)

#define OPEN_READ(x, y)                                 \
do {                                                    \
    x = fopen(y, "r");                                  \
    if(x == NULL)                                       \
    {                                                   \
        printf("%s failed to open for read\n", y);      \
        exit(-1);                                       \
    }                                                   \
} while(0)

#define OPEN_WRITE(x, y)                                \
do {                                                    \
    x = fopen(y, "w");                                  \
    if(x == NULL)                                       \
    {                                                   \
        printf("%s failed to open for write\n", y);     \
        exit(-1);                                       \
    }                                                   \
} while(0)

typedef double dtype;
//typedef MPI_DOUBLE MPI_DTYPE;

typedef struct my_matrix{
    dtype *storage;         /*storage for matrix data*/
    dtype **data;           /*matrix data*/
    int dims[NDIMS];        /*dimensions of full matrix*/
    int block_dims[NDIMS];  /*dimensions of current block, if sequential program block_dims == dims*/
} my_matrix;

/*Malloc that also checks if malloc failed*/
void *my_malloc_serial(int bytes);

/*Calloc that also checks if Calloc failed*/
void *my_calloc_serial(int num, int size);

/*Allocates space for matrix
    Args:   IN  - matrix_dims = {Number of rows, Number of columns}
            OUT - storage, matrix pointer to matrix storage
            OUT - A, pointer to matrix values
*/
void allocate_matrix(my_matrix *Matrix);

/*Allocates space for matrix will all values initialized to 0
    Args:   IN  - matrix_dims = {Number of rows, Number of columns}
            OUT - storage, matrix pointer to matrix storage
            OUT - A, pointer to matrix values
*/
void calloc_matrix(my_matrix *Matrix);

/*Frees the given my_matrix structure along with its data*/
void free_matrix(my_matrix **Matrix);

/*Writes the given matrix, A, into file_name
    Inputs: file_name: file to be written to
            Matrix: a my_matrix structure
 */
void write_matrix(char *file_name, my_matrix *Matrix);

/*Reads the graph from the given binary file name.
  Function will malloc space for the graph
    Args:   IN  - file_name: binary file to read from
            OUT - my_matrix structure, will be allocated in function
 */
void read_matrix(char *file_name, my_matrix *Matrix);

/*Will print matrix's num_node by num_node adjacency matrix
    Args:   IN  - preallocated my_matrix structure
 */
void print_matrix(my_matrix *Matrix);

#endif
