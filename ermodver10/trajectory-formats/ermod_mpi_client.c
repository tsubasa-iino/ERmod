#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>

#include "ermod_mpi_client.h"

static const char* server_conffile = "MpiConn";
static const char* default_service_name = "ermod-trajio";
static const int ERMOD_TAG_CELL = 1;
static const int ERMOD_TAG_COORD = 2;

typedef struct handler_t {
  MPI_Comm comm;
  int server_rank;
} handler;

static void show_mpi_error(int retval) {
  char buf[MPI_MAX_ERROR_STRING];
  int buflen;
  
  MPI_Error_string(retval, buf, &buflen);
  fprintf(stderr, "MPI returned unrecoverable error: \"%s\"\n", buf);
}

static void halt_with_mpi_error(int retval) {
  show_mpi_error(retval);
  exit(1);
}

static void sleep_real(double len) {
  int sec = (int)len;
  useconds_t usec = (len - sec) * 1000000;
  if(usec >= 1000000) usec = 1000000;
  sleep(sec);
  usleep(usec);
}

/*
  Connect to ERmod program via MPI intercommunicator.
 */
void* ermod_connect()
{
  char servicename[80];
  char port[MPI_MAX_PORT_NAME];
  handler* handle;
  int r;
  double timer = 0.1; /* 0.1 sec grace time */
  
  /* set servicename */
  do{
    FILE *fp;
    fp = fopen(server_conffile, "rt");
    if(fp == NULL) {
      /* failed to open => use default */
      strcpy(servicename, default_service_name);
      break;
    }
    fscanf(fp, "%s", servicename);
  }while(0);

  for(;;){
    /* Look up for MPI published name. */
    r = MPI_Lookup_name(servicename, MPI_INFO_NULL, port);
    if(r == 0) {
      /* Found port name */
      break;
    }
    if(r != MPI_ERR_NAME) {
      /* If name is not published, MPI_ERR_NAME is returned. 
	 Otherwise halt with error. */
      halt_with_mpi_error(r);
    }

    timer *= 1.2;
    if(timer > 10.0) timer = 10.0;
    sleep_real(timer); /* wait and retry */
  }

  handle = malloc(sizeof(handler));
  if(!handle) return NULL;

  /* 
     Note for future implementers: if you'd like to gather trajectory data,
     use MPI_COMM_WORLD instead of MPI_COMM_SELF and try to connect from all.
     You also have to modify the server.
  */
  r = MPI_Comm_connect(port, MPI_INFO_NULL, 0, MPI_COMM_SELF, &handle->comm);
  if(r){
    halt_with_mpi_error(r);
  }

  {
    int myrank, commsize;
    r = MPI_Comm_size(handle->comm, &commsize);
    if(r) halt_with_mpi_error(r);
    r = MPI_Comm_rank(handle->comm, &myrank);
    if(r) halt_with_mpi_error(r);
    if(commsize != 2) {
      fprintf(stderr, "Warning: ERmod MPI trajectory client: Communicator size is expected to be 2, but was %d.\n", commsize);
    }
    if(myrank != 1) {
      fprintf(stderr, "Warning: ERmod MPI trajectory client: Client rank is expected to be 1, but was %d.\n", myrank);
    }
    
    handle->server_rank = 1 - myrank;
  }

  return handle;
}

/*
  Close the connection to ERmod and release the resource.
 */
void ermod_disconnect(void* handle_p)
{
  handler* handle = (handler*)handle_p;

  MPI_Comm_disconnect(&handle->comm);
  MPI_Comm_free(&handle->comm);

  free(handle);
}

/*
  
 */
int ermod_send_trajctory(void* handle_p, int natoms, double cell[], double coords[])
{
  handler* handle = (handler*)handle_p;
  int r;
  /*
    Note: Turning this into MPI_Isend is a bit tricky, since it must have "clean-up" calls...
   */
  r = MPI_Send(&cell[0], 3 * 3, MPI_DOUBLE_PRECISION,
	       handle->server_rank, ERMOD_TAG_CELL,
	       handle->comm);
  if(r){ show_mpi_error(r); return 1; }
  r = MPI_Send(&coords[0], 3 * natoms, MPI_DOUBLE_PRECISION,
	       handle->server_rank, ERMOD_TAG_COORD,
	       handle->comm);
  if(r){ show_mpi_error(r); return 1; }
  return 0;
}

/* Fortran interfaces */

void ermod_connect_(void** handle_p_p)
{
  void* ptr = ermod_connect();
  *handle_p_p = ptr;
}

void ermod_disconnect_(void** handle_p_p)
{
  ermod_disconnect(*handle_p_p);
}

void ermod_send_trajctory_(void** handle_p_p, int* natoms_p, double cell[], double coords[], int* status_p)
{
  int status = ermod_send_trajctory(handle_p_p, *natoms_p, cell, coords);
  *status_p = status;
}

