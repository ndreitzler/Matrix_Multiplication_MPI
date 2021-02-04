//Written by Dr. Walt Ligon of Clemson University 

// #include <stdio.h>
// #include <stdlib.h>
//#include <mpi.h>
#include "graph_checkerboard_io.h"
//#include "matrix.h"

void error_out(int ret, int ID, MPI_Status *status)
{
    if (ret != MPI_SUCCESS)
    {
        char string[MPI_MAX_ERROR_STRING];
        int strsz;
        MPI_Error_string(ret, string, &strsz);
        debug("%d: error: %s\n", ID );
        MPI_Finalize();
        exit(-1);
    }
    ret = MPI_SUCCESS;
}

/*
 *   Function 'read_checkerboard_graph' reads a graph from
 *   a file. The first element of the file is an integer
 *   whose value is the number of nodes in the graph, and
 *   the size each size of the adjacencymatrix ('n' rows
 *   and 'n' columns). What follows are 'n'*'n' values
 *   representing the matrix elements stored in row-major
 *   order.  This function allocates blocks of the matrix to
 *   the MPI processes.
 *
 *   The number of processes must be a square number.
 *   You must pass in communicator with a 2D grid topology.
 *   This routines allocates and frees the 2D matrix partition
 */
 
void read_checkerboard_graph (
                char *s,              /* IN - File name */
                my_matrix *Matrix,
                // void ***subs,         /* OUT - 2D array */
                // void *Matrix->storage,       /* OUT - Array elements */
                MPI_Datatype dtype_this,   /* IN - Element type */
                // int *N,               /* OUT - Array cols */
                MPI_Comm grid_comm)   /* IN - Communicator */
{
    //void *** subs = Matrix->
    int          datum_size = 0;      /* Bytes per elements */
    int          grid_coord[2] = {0, 0};  /* Process coords */
    int          grid_id = 0;         /* Process rank */
    int          grid_period[2] = {0, 0}; /* Wraparound */
    int          grid_size[2] = {0, 0};   /* Dimensions of grid */
    int          i, ID = 0, ret = MPI_SUCCESS;
    void       **lptr = NULL;         /* Pointer into 'subs' */
    void        *rptr = NULL;         /* Pointer into 'storage' */
    int          matsize[2] = {0, 0};
    int          subsizes[2] = {0, 0};
    int          starts[2] = {0, 0};
    MPI_Status   status;              /* Results of read */
    MPI_File     fh;
    MPI_Datatype subarray_t;

    MPI_Comm_rank (grid_comm, &grid_id); ID = grid_id;

    debug( "%d: inside read_checkboard_graph\n", ID );

    MPI_Type_size(dtype_this, &datum_size);

    debug( "%d: grid_id = %d dtype_this_size = %d\n", ID, grid_id, datum_size );

    ret = MPI_Cart_get (grid_comm, 2, grid_size, grid_period, grid_coord);
    error_out(ret, ID, NULL);

    debug( "%d: g_size[0] = %d, g_size[1] = %d\n", ID, grid_size[0], grid_size[1] );
    debug( "%d: g_peri[0] = %d, g_peri[1] = %d\n", ID, grid_period[0], grid_period[1] );
    debug( "%d: g_coor[0] = %d, g_coor[1] = %d\n", ID, grid_coord[0], grid_coord[1] );

    ret = MPI_File_open(grid_comm, s, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    debug( "%d: Open file: %s\n", ID, s );
    error_out(ret, ID, NULL);

    if (grid_id == 0)
    {
        //Read first 2 elements
        //ret = MPI_File_read(fh, N, 1, MPI_INTEGER, &status);
        ret = MPI_File_read(fh, Matrix->dims, 2, MPI_INTEGER, &status);
        error_out(ret, ID, &status);
        debug( "%d: Read size from file: %s\n", ID, s );
    }
    MPI_Bcast(Matrix->dims, 2, MPI_INTEGER, 0, grid_comm);
    debug( "%d: N = %d\n", ID, *N );

    matsize[0] = Matrix->dims[0];
    matsize[1] = Matrix->dims[1];
    debug( "%d: matsize[0] = %d, matsize[1] = %d\n", ID, matsize[0], matsize[1] );

    subsizes[0] = BLOCK_SIZE(grid_coord[0], grid_size[0], matsize[0]);
    subsizes[1] = BLOCK_SIZE(grid_coord[1], grid_size[1], matsize[1]);
    debug( "%d: subsz[0] = %d, subsz[1] = %d\n", ID, subsizes[0], subsizes[1] );

    Matrix->block_dims[0] = subsizes[0];
    Matrix->block_dims[1] = subsizes[1];

    starts[0] = BLOCK_LOW(grid_coord[0], grid_size[0], matsize[0]);
    starts[1] = BLOCK_LOW(grid_coord[1], grid_size[1], matsize[1]);
    debug( "%d: starts[0] = %d, starts[1] = %d\n", ID, starts[0], starts[1] );

    /* Dynamically allocate two-dimensional matrix 'subs' */

    allocate_matrix(Matrix);

    // Matrix->storage = (void *) malloc (subsizes[0] * subsizes[1] * datum_size);
    // Matrix->data = (void **) malloc (subsizes[0] * sizeof(void *));
    // lptr = (void **) Matrix->data;
    // rptr = (void *) Matrix->storage;
    // for (i = 0; i < subsizes[0]; i++)
    // {
    //     *(lptr++) = (void *) rptr;
    //     rptr += (subsizes[1] * datum_size);
    // }

    debug( "%d: create_subarray\n", ID );
    MPI_Type_create_subarray(2, matsize, subsizes, starts, MPI_ORDER_C, dtype_this, &subarray_t);

    debug( "%d: commit\n", ID );
    MPI_Type_commit(&subarray_t);

    debug( "%d: set_view\n", ID );
    MPI_File_set_view(fh, sizeof(dtype_this), dtype_this, subarray_t, "native", MPI_INFO_NULL);

    debug( "%d: read file data\n", ID );
    ret = MPI_File_read(fh, Matrix->storage, subsizes[0] * subsizes[1], dtype_this, &status);
    error_out(ret, ID, &status);

    debug( "%d: close\n", ID );
    MPI_File_close(&fh);

    debug( "%d: free\n", ID );
    MPI_Type_free(&subarray_t);

    debug( "%d: read returning\n", ID );
}

/*
 *   Function 'write_checkerboard_graph' writes a graph to
 *   a file. The first element of the file is an integer
 *   whose value is the number of nodes in the graph, and
 *   the size each size of the adjacency matrix ('n' rows
 *   and 'n' columns). What follows are n x n values
 *   representing the matrix elements stored in row-major
 *   order.  This function allocates frees the matrix for
 *   the MPI processes.
 *
 *   The number of processes must be a square number.
 */
 
void write_checkerboard_graph (
                char *s,              /* IN - File name */
                //void ***subs,         /* OUT - 2D array */
                //void *Matrix->storage,       /* OUT - Array elements */
                my_matrix *Matrix,
                MPI_Datatype dtype_this,   /* IN - Element type */
                //int *N,                /* OUT - Array cols */
                MPI_Comm grid_comm)   /* IN - Communicator */
{
    int          datum_size = 0;      /* Bytes per elements */
    int          grid_coord[2] = {0, 0};  /* Process coords */
    int          grid_id = 0;         /* Process rank */
    int          grid_period[2] = {0, 0}; /* Wraparound */
    int          grid_size[2] = {0, 0};   /* Dimensions of grid */
    int          i, ID = 0, ret = 0;
    void       **lptr = NULL;         /* Pointer into 'subs' */
    void        *rptr = NULL;         /* Pointer into 'storage' */
    int          matsize[2] = {0, 0};
    int          subsizes[2] = {0, 0};
    int          starts[2] = {0, 0};
    MPI_Status   status;              /* Results of read */
    MPI_File     fh;
    MPI_Datatype subarray_t;

    MPI_Comm_rank (grid_comm, &grid_id); ID = grid_id;

    debug( "%d: inside write_checkboard_graph\n", ID );

    MPI_Type_size(dtype_this, &datum_size);

    debug( "%d: grid_id = %d dtype_this_size = %d\n", ID, grid_id, datum_size );

    ret = MPI_Cart_get (grid_comm, 2, grid_size, grid_period, grid_coord);
    if (ret != MPI_SUCCESS)
    {
        char string[MPI_MAX_ERROR_STRING];
        int strsz;
        MPI_Error_string(ret, string, &strsz);
        //debug("%d: error: %s\n", ID );
        MPI_Finalize();
        exit(-1);
    }

    debug( "%d: g_size[0] = %d, g_size[1] = %d\n", ID, grid_size[0], grid_size[1] );
    debug( "%d: g_peri[0] = %d, g_peri[1] = %d\n", ID, grid_period[0], grid_period[1] );
    debug( "%d: g_coor[0] = %d, g_coor[1] = %d\n", ID, grid_coord[0], grid_coord[1] );
    debug( "%d: *N = %d\n", ID, *N );

    ret = MPI_File_open(grid_comm, s, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    error_out(ret, ID, NULL);
    debug( "%d: Open file: %s\n", ID, s );

    if (grid_id == 0)
    {
        //printf("%d: %d %d\n", grid_id, N[0], N[1]);
        ret = MPI_File_write(fh, Matrix->dims, 2, MPI_INTEGER, &status);
        //fwrite(N, sizeof(int), 2, fh);
        //ret = MPI_File_write(fh, &M, 1, MPI_INTEGER, &status);
        error_out(ret, ID, &status);
        debug( "%d: write size to file: %s\n", ID, s );
    }

    matsize[0] = Matrix->dims[0];
    matsize[1] = Matrix->dims[1];
    debug( "%d: matsize[0] = %d, matsize[1] = %d\n", ID, matsize[0], matsize[1] );

    subsizes[0] = BLOCK_SIZE(grid_coord[0], grid_size[0], matsize[0]);
    subsizes[1] = BLOCK_SIZE(grid_coord[1], grid_size[1], matsize[1]);
    debug( "%d: subsz[0] = %d, subsz[1] = %d\n", ID, subsizes[0], subsizes[1] );

    starts[0] = BLOCK_LOW(grid_coord[0], grid_size[0], matsize[0]);
    starts[1] = BLOCK_LOW(grid_coord[1], grid_size[1], matsize[1]);
    debug( "%d: starts[0] = %d, starts[1] = %d\n", ID, starts[0], starts[1] );

    debug( "%d: create_subarray\n", ID );
    MPI_Type_create_subarray(2, matsize, subsizes, starts, MPI_ORDER_C, dtype_this, &subarray_t);

    debug( "%d: commit\n", ID );
    MPI_Type_commit(&subarray_t);
    
    debug( "%d: set_view\n", ID );
    MPI_File_set_view(fh, sizeof(dtype_this), dtype_this, subarray_t, "native", MPI_INFO_NULL);

    debug( "%d: write data to file\n", ID );
    ret = MPI_File_write(fh, Matrix->storage, subsizes[0] * subsizes[1], dtype_this, &status);
    error_out(ret, ID, &status);

    debug( "%d: close file\n", ID );
    ret = MPI_File_close(&fh);
    error_out(ret, ID, NULL);

    debug( "%d: free subtype\n", ID );
    ret = MPI_Type_free(&subarray_t);
    error_out(ret, ID, NULL);

    // free(Matrix->storage);
    // Matrix->storage = NULL;
    // free(*subs);
    // *subs = NULL;

    debug( "%d: write returning\n", ID );
}
