/*  Name: Nick Dreitzler
    Username: ndreitz
    CUID: C15564443
    Class: ECE 4730
    Professor: Walt Ligon
    Semester: Fall2020
    Filename: floyd-parallel.c

    Purpose: Reads 2 matrix files. Multiplies the first matrix
        the second matrix. Writes resulting matrix to a file.
        Assumes all matrices are square.

    usage: ./mm-parallel <input_file1>  <input_file2> <output_file>
    
    Output: Binary file with n by n matrix representing 
        Floyd's shortest pairs.
 */

#include "matrix.h"
#include "graph_checkerboard_io.h"
//#include "MyMPI.h"

#define NDIMS 2
#define R_ROW 0
#define R_COL 1
#define R_WORLD 0

const MPI_Datatype MPI_DTYPE = MPI_DOUBLE;

void matrix_multiply_parallel( my_matrix *Matrix1, my_matrix *Matrix2, my_matrix *Result, MPI_Comm Grid_comm);
void copy_storage(dtype *source_storage, dtype *copier_storage, int size);
void matrix_mult_part( my_matrix *Matrix1, my_matrix *Matrix2, my_matrix *Result);
int find_sqrt(int size);

int main(int argc, char **argv)
{
    MPI_Comm MPI_grid;  /*MPI communicator with cartesian topology*/
    my_matrix *Matrix1 = (my_matrix *)malloc(sizeof(my_matrix));
    my_matrix *Matrix2 = (my_matrix *)malloc(sizeof(my_matrix));
    my_matrix *Result = (my_matrix *)malloc(sizeof(my_matrix));
    double exe_time;    /*Execution time without IO*/
    double Max_exe_time; /*Max execution time without IO*/
    //double TimeIO;      /*Execution time with IO*/
    //double MaxTimeIO;   /*Max execution time with IO*/
    int size;           /*Number of MPI tasks*/
    int rank;           /*rank of MPI task*/
    int dims[] = {0,0}; /*Dimensions of MPI_grid*/
    int periods[NDIMS]; /*Periods of MPI_grid*/

    //MPI setup
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Need input and output files
    if(argc != 4)
    {
        if(rank == 0)
            printf("Usage:\n\t./floyd-serial <input file1> <input_file2> <output file>\n");
        MPI_Finalize();
        exit(1);
    }

    //Create cartesian topology and read adjacency matrix into tasks*/
    MPI_Dims_create( size , NDIMS , dims);
    MPI_Cart_create( MPI_COMM_WORLD , NDIMS , dims , periods , 1 , &MPI_grid);

    //Read in both matrices
    read_checkerboard_graph(argv[1], Matrix1, MPI_DTYPE, MPI_grid);
    read_checkerboard_graph(argv[2], Matrix2, MPI_DTYPE, MPI_grid);

    exe_time = -MPI_Wtime();
    matrix_multiply_parallel(Matrix1, Matrix2, Result, MPI_grid);
    exe_time += MPI_Wtime();

    //Write resulting matrix to file
    write_checkerboard_graph(argv[3], Result, MPI_DTYPE, MPI_grid);

    //Determine max execution time and store result with task 0
    //MPI_Reduce( (void *)(&TimeIO), (void *)(&MaxTimeIO), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce( (void *)(&exe_time), (void *)(&Max_exe_time), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Print results
    // if(rank == 0)
    // {
    //     printf("%5s\t%5s\t%11s\t%13s\n", "size", "n", "TimeIO (ms)", "exe_time (ms)");
    //     printf("%5d\t%5d\t%11f\t%13f\n", size, num_nodes, TimeIO*1000, exe_time*1000);
    // }

    // Print results
    if(rank == 0)
    {
        printf("%5s\t%5s\t%5s\t%13s\n", "","size", "n", "exe_time (ms)");
        printf("%5s\t%5d\t%5d\t%13f\n", "data",size, Result->dims[ROW], exe_time*1000);
    }

    //Clean up
    free_matrix(&Result);
    free_matrix(&Matrix1);
    free_matrix(&Matrix2);
    MPI_Comm_free( &MPI_grid);
    MPI_Finalize();
    return(0);
}

/*Multiplies matrix 1 by matrix 2.
    Assumes that there is a square number of processors.
*/
void matrix_multiply_parallel(
    my_matrix *Matrix1,     /*IN  - multiplicand*/
    my_matrix *Matrix2,     /*IN  - multiplier*/
    my_matrix *Result,      /*OUT - result of multiplication*/
    MPI_Comm Grid_comm)     /*IN  - MPI communicator with cartesian topology*/
{
    MPI_Status status;      /*MPI status*/
    my_matrix *tmp_block = (my_matrix *)malloc(sizeof(my_matrix)); /*tmp block used for message passing*/
    int storage_size;       /*Number of dtype in a block*/
    int dims[NDIMS];        /*Holds dimensions of Grid_comm*/
    int periods[NDIMS];     /*Holds periods of Grid_comm*/
    int grid_coords[NDIMS]; /*Holds coordnates of node*/
    int block_size[NDIMS];  /*dimensions of current block, {number of rows, number of columns}*/
    int send_rank;          /*rank to send data too*/
    int recv_rank;          /*rank to recieve data from*/
    int rank;               /*processor rank in MPI_COMM_WORLD*/
    int size;               /*number of procs in MPI_COMM_WORLD, must be a perfect square*/
    int sqrt_size;          /*sqrt of size, must be an int*/
    int i;

    //Get MPI data
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Cart_get( Grid_comm , NDIMS , dims , periods , grid_coords);

    sqrt_size = find_sqrt(size);

    //Check that size is a perfect square
    if(sqrt_size == -1)
    {
        if(rank == 0)
            printf("Number of procs must be a perfect square\n");
        MPI_Abort (MPI_COMM_WORLD, -1);
    }

    //Check that number of rows is divisable by sqrt(size)
    if(!(Matrix1->dims[ROW] % sqrt_size == 0))
    {
        if(rank == 0)
            printf("Elements in a dimension must be divisable by sqrt of size\n");
        MPI_Abort (MPI_COMM_WORLD, -1);
    }

    //Check dimensions match
    if(Matrix1->dims[COL] != Matrix2->dims[ROW])
    {
        if(rank == 0)
            printf("matrix 1 columns must match matrix 2 rows\n");
        MPI_Abort (MPI_COMM_WORLD, -1);
    }

    //Check matruix 1 is square
    if(Matrix1->dims[ROW] != Matrix1->dims[COL])
    {
        if(rank == 0)
            printf("matrix 1 must be a square matrix\n");
        MPI_Abort (MPI_COMM_WORLD, -1);
    }

    //Check matrix 2 is square
    if(Matrix2->dims[ROW] != Matrix2->dims[COL])
    {
        if(rank == 0)
            printf("matrix 2 cmust be a square matrix\n");
        MPI_Abort (MPI_COMM_WORLD, -1);
    }

    //Set matrix dims for full result. not just current block
    Result->dims[ROW] = Matrix1->dims[ROW];
    Result->dims[COL] = Matrix2->dims[COL];

    //Determine number of rows/cols that the the current task has.
    //Program assumes that number of elements is a row/col is divisible
    //by the number of processors.
    block_size[ROW] = BLOCK_SIZE(grid_coords[ROW], dims[ROW], Result->dims[ROW]);
    block_size[COL] = BLOCK_SIZE(grid_coords[COL], dims[COL], Result->dims[COL]);

    tmp_block->block_dims[ROW] = block_size[ROW];
    tmp_block->block_dims[COL] = block_size[COL];

    Result->block_dims[ROW] = block_size[ROW];
    Result->block_dims[COL] = block_size[COL];

    storage_size = tmp_block->block_dims[ROW]*tmp_block->block_dims[COL];

    //Allcate space for tmp block used for sending
    allocate_matrix(tmp_block);

    //Calloc space for results block.
    //Use calloc so values are initalized to 0
    calloc_matrix(Result);

    //Align phase
    //Procs not in the first row
    if(grid_coords[ROW])
    {
        //Copy data into tmp_block then send_recv
        copy_storage(Matrix1->storage, tmp_block->storage, storage_size);
        //Shift matrix 1 blocks left once per row value
        MPI_Cart_shift( Grid_comm , 1 , -grid_coords[ROW] , &recv_rank , &send_rank);
        MPI_Sendrecv( (void *)(tmp_block->storage), storage_size , MPI_DTYPE , send_rank , 0 , (void *)(Matrix1->storage) , storage_size , MPI_DTYPE , recv_rank , 0 , MPI_COMM_WORLD , &status);
    }
    //Procs not in first column
    if(grid_coords[COL])
    {
        copy_storage(Matrix2->storage, tmp_block->storage, storage_size);
        //Shift matrix 2 blocks down once per column value
        MPI_Cart_shift( Grid_comm , 0 , -grid_coords[COL] , &recv_rank , &send_rank);
        MPI_Sendrecv( (void *)(tmp_block->storage), storage_size , MPI_DTYPE , send_rank , 0 , (void *)(Matrix2->storage) , storage_size , MPI_DTYPE , recv_rank , 0 , MPI_COMM_WORLD , &status);
    }

    //First matrix multiply
    matrix_mult_part(Matrix1, Matrix2, Result);

    for(i = 0; i < sqrt_size-1; i++)
    {
        //Shift matrix 1 left one column
        copy_storage(Matrix1->storage, tmp_block->storage, storage_size);
        MPI_Cart_shift( Grid_comm , 1 , -1 , &recv_rank , &send_rank);
        MPI_Sendrecv( (void *)tmp_block->storage, storage_size , MPI_DTYPE , send_rank , 0 , (void *)(Matrix1->storage) , storage_size , MPI_DTYPE , recv_rank , 0 , MPI_COMM_WORLD , &status);

        //Shift matrix 2 down one row
        copy_storage(Matrix2->storage, tmp_block->storage, storage_size);
        MPI_Cart_shift( Grid_comm , 0 , -1 , &recv_rank , &send_rank);
        MPI_Sendrecv( (void *)tmp_block->storage, storage_size , MPI_DTYPE , send_rank , 0 , (void *)(Matrix2->storage) , storage_size , MPI_DTYPE , recv_rank , 0 , MPI_COMM_WORLD , &status);
    
        //Subsequent matrix multiple. sqrt_size in total
        matrix_mult_part(Matrix1, Matrix2, Result);
    }//End for i

    //Clean up
    free_matrix(&tmp_block);
}

/*Copies the data in source_storage to copier_storage*/
void copy_storage(dtype *source_storage, dtype *copier_storage, int size)
{
    int i;
    for(i = 0; i < size; i++)
            copier_storage[i] = source_storage[i];
}

/*Multiplies two square matrices together and stores the result in a new matrix.
    Assumes that the result matrix exists and has been allocated space
*/
void matrix_mult_part(
    my_matrix *Matrix1,     /*IN  - multiplicand*/
    my_matrix *Matrix2,     /*IN  - multiplier*/
    my_matrix *Result)      /*IN/OUT - result of multiplication*/    
{
    int i,j,k;

    for (i = 0; i < Result->block_dims[ROW]; ++i)
        for (j = 0; j < Result->block_dims[ROW]; ++j)
            for (k = 0; k < Result->block_dims[ROW]; ++k)
                Result->data[i][j] += Matrix1->data[i][k] * Matrix2->data[k][j];
}


/*Finds the sqrt of size if size is between 1-64 and a perfect square*/
int find_sqrt(int size)
{
    int out = -1;

    switch (size)
    {
    case 1:
        out = 1;
        break;
    case 4:
        out = 2;
        break;
    case 9:
        out = 3;
        break;
    case 16:
        out = 4;
        break;
    case 25:
        out = 5;
        break;
    case 36:
        out = 6;
        break;
    case 49:
        out = 7;
        break;
    case 64:
        out = 8;
        break;
    default:
        break;
    }
    return out;
}

// //Send current b vector to tasks in first row
//         root = BLOCK_OWNER(k, dims[COL], Result->dims[COL]);
//         //if(ranks[R_WORLD] == 0)
//         //    printf("%d: %d\n", ranks[R_WORLD], root);
//         if(root == ranks[R_COL])
//         {
//             //Task at {0,0} would send to itself
//             if((ranks[R_ROW] != 0) || (ranks[R_COL] != 0))
//             {
//                 for(i = 0; i < tmp_block->dims[ROW]; i++)
//                 {
//                     //BLOCK_LOW(root, dims[1], num_nodes);
//                     offset = k - BLOCK_LOW(root, dims[COL], Result->dims[ROW]);
//                     printf("%d: %d \n", ranks[R_COL], offset);
//                     tmp_b_send[i] = M2[i][offset];
//                 }
//                 target_coords[0] = 0;
//                 target_coords[1] = ranks[R_ROW];
//                 MPI_Cart_rank( Grid_comm , target_coords , &target_rank);
//                 MPI_Send( tmp_b_send , tmp_block->dims[ROW] , MPI_DTYPE, target_rank , 0 , Grid_comm);
//             } else {
//                 for(i = 0; i < tmp_block->dims[ROW]; i++)
//                 {
//                     tmp_b[i] = M2[i][k];
//                 }
//             }
//         }//End if(root == ranks[R_COL])
        
//         if((ranks[R_ROW] == 0) && ranks[R_COL] != 0)
//         {
//             target_coords[0] = ranks[R_COL];
//             target_coords[1] = root;
//             MPI_Cart_rank( Grid_comm , target_coords , &target_rank);
//             MPI_Recv( tmp_b , tmp_block->dims[COL] , MPI_DTYPE, target_rank , 0 , Grid_comm, &status);

//             sleep(ranks[0]);
//             for(i = 0; i < tmp_block->dims[R_COL]; i++)
//             {
//                 printf("%0.6f ", tmp_b[i]);
//                 ffluy("\n");
//             fflush(stdout);
//         }



        //         if(grid_coords[ROW] == 1)
        // {
        //     sleep(ranks[R_WORLD]);
        //     for(i = 0; i < tmp_block->dims[ROW]; i++)
        //     {
        //         for(j = 0; j < tmp_block->dims[COL]; j++)
        //                 printf("%0.6f ", matrix1[i][j]);
        //         printf("---%d\n",grid_coords[COL]);
        //     }

        // }