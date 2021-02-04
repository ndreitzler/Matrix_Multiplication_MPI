#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "matrix.h"


/*
 *   Function 'read_checkerboard_graph' reads a graph from
 *   a file. The first element of the file is an integer
 *   whose value is the number of nodes in the graph, and
 *   the size of each side of the adjacency matrix ('n' rows
 *   and 'n' columns). What follows are 'n'*'n' integer values
 *   representing the matrix elements stored in row-major
 *   order.  This function allocates blocks of the matrix to
 *   the MPI processes, the caller need only pass in the pinters.
 *
 *   Function write_checkboard_graph does the exact opposite,
 *   writing (overwriting or creating) a file with the matrix as
 *   described above and then frees the buffers previously allocated.
 *
 *   Output oarguments of each function must be declared as shown.
 *   The MPI Communicator must have a 2D Cartesian topology.
 *
 *   The number of processes must be a square number.
 *   (This can actually be relaxed to a generic 2D grid)
 *
 *   Defines in this file can be used to turn on extensive debugging
 *   messages, and can be used for message in the calling code.
 *   Block macros are repeated here for convenience.
 */

//typedef MPI_DOUBLE MPI_DTYPE
//const MPI_Datatype MPI_DTYPE_MM = MPI_DOUBLE;

#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p, n) -1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id)+1, p, n)-BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p)*((index)+1)-1)/(n))

int get_size (MPI_Datatype t);

void *my_malloc (
   int id,     /* IN - Process rank */
   int bytes);  /* IN - Bytes to allocate */

void read_checkerboard_graph (
   char *s,              /* IN - File name */
   my_matrix *Matrix,
   // void ***subs,         /* OUT - 2D array */
   // void *Matrix->storage,       /* OUT - Array elements */
   MPI_Datatype dtype_this,   /* IN - Element type */
   // int *N,               /* OUT - Array cols */
   MPI_Comm grid_comm);   /* IN - Communicator */

void write_checkerboard_graph (
   char *s,              /* IN - File name */
   //void ***subs,         /* OUT - 2D array */
   //void *Matrix->storage,       /* OUT - Array elements */
   my_matrix *Matrix,
   MPI_Datatype dtype_this,   /* IN - Element type */
   //int *N,                /* OUT - Array cols */
   MPI_Comm grid_comm);   /* IN - Communicator */

void error_out(int ret, int ID, MPI_Status *status);

#define NOT_DEBUG
#ifdef DEBUG
#define debug(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)
#else
#define debug(fmt, ...)
#endif

