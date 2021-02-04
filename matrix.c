/*  Name: Nick Dreitzler
    Username: ndreitz
    CUID: C15564443
    Class: ECE 4730
    Professor: Walt Ligon
    Semester: Fall2020
    Filename: make-matrix.c

    Purpose: Various functions to read/write/print matrix
 */

#include "matrix.h"

/*Malloc that also checks if malloc failed*/
void *my_malloc_serial (int bytes)
{
   void *buffer;
   if ((buffer = malloc ((size_t) bytes)) == NULL) {
      printf ("Error: Malloc failed\n");
      exit(1);
   }
   return buffer;
}

/*Malloc that also checks if malloc failed*/
void *my_calloc_serial (int num, int size)
{
   void *buffer;
   if ((buffer = calloc (num, size)) == NULL) {
      printf ("Error: Calloc failed\n");
      exit(1);
   }
   return buffer;
}

/*Allocates space for matrix adjacency matrix
    Args:   IN  - Matrix->dims = {Number of rows, Number of columns}
            OUT - S, matrix pointer to matrix storage
            OUT - A, pointer to matrix values
*/
void allocate_matrix(my_matrix *Matrix)
{
    int i;

    //printf("%d %d\n", Matrix->block_dims[ROW], Matrix->block_dims[COL]);

    //Allocate space for matrix
    Matrix->storage = (dtype *)my_malloc_serial(sizeof(dtype)*Matrix->block_dims[ROW]*Matrix->block_dims[COL]);
    Matrix->data = (dtype **)my_malloc_serial(sizeof(dtype *)*Matrix->block_dims[ROW]);
    for(i = 0; i < Matrix->block_dims[ROW]; i++)
        Matrix->data[i] = &(Matrix->storage[Matrix->block_dims[COL]*i]);

}

/*Allocates space for matrix adjacency matrix using calloc for storage
    Args:   IN  - Matrix->dims = {Number of rows, Number of columns}
            OUT - S, matrix pointer to matrix storage
            OUT - A, pointer to matrix values
*/
void calloc_matrix(my_matrix *Matrix)
{
    int i;

    //Allocate space for matrix
    Matrix->storage = (dtype *)my_calloc_serial(sizeof(dtype),Matrix->block_dims[ROW]*Matrix->block_dims[COL]);
    Matrix->data = (dtype **)my_malloc_serial(sizeof(dtype *)*Matrix->block_dims[ROW]);
    for(i = 0; i < Matrix->block_dims[ROW]; i++)
        Matrix->data[i] = &(Matrix->storage[i*Matrix->block_dims[COL]]);
}

/*Frees the given my_matrix structure along with its data*/
void free_matrix(my_matrix **Matrix)
{
    //printf("test1");
    FREE((*Matrix)->storage);
    FREE((*Matrix)->data);
    FREE(*Matrix);
}


/*Writes the given matrix, A, into file_name
    Inputs: file_name: file to be written to
            Matrix: a my_matrix structure
 */
void write_matrix(char *file_name, my_matrix *Matrix)
{
    FILE *fptr;
    int i;
    
    OPEN_WRITE(fptr, file_name);

    //Write number of rows and number of columns to file
    fwrite(Matrix->dims, 2, sizeof(int), fptr);

    //write 1 row at a time
    for (i = 0; i < Matrix->dims[ROW]; i++)
    {
        fwrite(Matrix->data[i], Matrix->dims[COL], sizeof(dtype), fptr);
    }

    CLOSE(fptr);
}

/*Reads the graph from the given binary file name.
  Function will malloc space for the graph
    Args:   IN  - file_name: binary file to read from
            OUT - my_matrix structure, will be allocated in function
 */
void read_matrix(char *file_name, my_matrix *Matrix)
{
    FILE *fptr;
    // dtype *storage;
    // dtype **matrix;
    int i;

    OPEN_READ(fptr, file_name);

    //Read number of rows and number of columns
    fread(Matrix->dims, 2, sizeof(int), fptr);

    Matrix->block_dims[ROW] = Matrix->dims[ROW];
    Matrix->block_dims[COL] = Matrix->dims[COL];

    allocate_matrix(Matrix);

    //Read contents
    for (i = 0; i < Matrix->dims[ROW]; i++)
    {
        fread(Matrix->data[i], Matrix->dims[COL], sizeof(dtype), fptr);
    }

    CLOSE(fptr);

    // *S = storage;
    // *A = matrix;
}

/*Will print matrix's num_node by num_node adjacency matrix
    Args:   IN  - preallocated my_matrix structure
 */
void print_matrix(my_matrix *Matrix)
{
    int i, j;

    // Line 1
    printf("%4s", "|");
    for(i = 0; i < Matrix->block_dims[COL]; i++)
        printf("%14d", i);
    printf("\n");

    //Line 2
    printf("%4s", "|");
    for(i = 0; i < Matrix->block_dims[COL]; i++)
        printf("%s", "--------------");
    printf("-\n");

    //The rest
    for(i = 0; i < Matrix->block_dims[ROW]; i++)
    {
        printf("%2d%2s", i, "|");
        for(j = 0; j < Matrix->block_dims[COL]; j++)
        {
            printf("%14.6f", Matrix->data[i][j]);
        }
        printf("\n");
    }
}

