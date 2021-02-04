/*  Name: Nick Dreitzler
    Username: ndreitz
    CUID: C15564443
    Class: ECE 4730
    Professor: Walt Ligon
    Semester: Fall2020
    Filename: make-matrix.c

    Purpose: To create a binary file of a n by n matrix. 
        

    Input:  -r <rows> number of rows in the matrix, default 100
            -c <cols> number of columns in the matrix, default 100
            -l <low> lower bound on the values in the matrix (lower case ell), default 50
            -u <up> upper bound of the values in the matrix, default 500000
            -o <output file>, default matrix.bin
            -h <print usage>
    
    Output: Binary file with n by n matrix
 */

#include "matrix.h"

#define BUFFER_SIZE 4096

void create_matrix(int *bounds, int *dims, my_matrix *Matrix);

int main(int argc, char **argv)
{
    my_matrix *Matrix = (my_matrix *)malloc(sizeof(my_matrix));
    // dtype *storage;         /*pointer all data for matrix*/
    // dtype **matrix;          /*Data written to file*/ 
    int dims[] = {100, 100};
    int bounds[] = {50, 500000}; /*{Lower bound of values in matrix, Upper bound}*/
    int opt;
    char filename_default[] = "matrix.bin";
    char *out_file = NULL;  /*Output file*/

    srand(time(0));

    //Parse Input
    while((opt = getopt(argc, argv, "r:c:l:u:o:h")) != -1) 
    {
        switch(opt) 
        { 
            case 'r':
                dims[ROW] = atoi(optarg);
                break;
            case 'c':
                dims[COL] = atoi(optarg);
                break;
            case 'l':
                bounds[LOWER] = atoi(optarg);
                break;
            case 'u':
                bounds[UPPER] = atoi(optarg);
                break;
            case 'o':
                out_file = strdup(optarg);
                break;
            case 'h':
                printf("Usage:\n");
                printf("\t-r <rows> number of rows in the matrix, default 100\n");
                printf("\t-c <cols> number of columns in the matrix, default 100\n");
                printf("\t-l <low> lower bound on the values in the matrix (lower case ell), default 50\n");
                printf("\t-u <up> upper bound of the values in the matrix, default 500000\n");
                printf("\t-o <output file>, default matrix.bin\n");
                printf("\t-h <print usage>\n");
                exit(1);
                break;
        }
    }

    //If no output file given use default
    if(out_file == NULL)
        out_file = strdup(filename_default);
        
    //Generate matrix
    create_matrix(bounds, dims, Matrix);

    //Will free matrix but not storage
    write_matrix(out_file, Matrix);

    //clean up
    free_matrix(&Matrix);
    FREE(out_file);

    return(0);
}

/*Creates a Matrix->dims[ROW] by Matrix->dims[COL] matrix. 
    Will allocate memory for the matrix.
  Args:     IN  - bounds = {lower bound, upper bound}
            IN  - matrix dims = {number of rows, number of columns}
            OUT - Matrix, my_matrix structrue
 */ 
void create_matrix(int *bounds, int *dims, my_matrix *Matrix)
{
    dtype temp; 
    int i,j;

    Matrix->dims[ROW] = dims[ROW];
    Matrix->dims[COL] = dims[COL];

    Matrix->block_dims[ROW] = dims[ROW];
    Matrix->block_dims[COL] = dims[COL];

    allocate_matrix(Matrix);

    for(i = 0; i < Matrix->dims[ROW]; i++)
    {
        for(j = 0; j < Matrix->dims[COL]; j++)
        {
            temp = ((double)bounds[LOWER] + ((double)rand()/(bounds[UPPER]-bounds[LOWER])))/1000.0;
            Matrix->data[i][j] = temp;
        }//End for j
    }//End for i
}
