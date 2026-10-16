/* ------- file: -------------------------- parallel.c --------------
       --------------------------                      ----------RH-- */

/* --- Splits a loop over N independent items into contiguous chunks
       and runs them concurrently, one chunk per thread. Meant for
       loops over grid points, lines, etc. whose iterations do not
       share state; the chunks are static, so this suits loops whose
       items cost about the same. The calling thread does chunk 0
       itself. With nthreads <= 1 (or N <= 1) it simply calls work()
       once for the whole range, with no threads involved.
       --                                              -------------- */

#include <stdlib.h>
#include <pthread.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"

extern InputData input;
extern char messageStr[];

typedef struct {
  ParallelRange work;
  void         *arg;
  int           begin, end, nt;
} Chunk;

static void *runChunk(void *argument)
{
  Chunk *c = (Chunk *) argument;

  c->work(c->begin, c->end, c->nt, c->arg);
  return NULL;
}

void parallelFor(int N, int nthreads, ParallelRange work, void *arg)
{
  int       nt, size, Nchunk;
  Chunk    *chunk;
  pthread_t *id;
  bool_t   *started;

  if (nthreads > N) nthreads = N;
  if (nthreads <= 1) {
    work(0, N, 0, arg);
    return;
  }

  size   = (N + nthreads - 1) / nthreads;
  Nchunk = (N + size - 1) / size;
  chunk  = (Chunk *) malloc(Nchunk * sizeof(Chunk));
  id     = (pthread_t *) malloc(Nchunk * sizeof(pthread_t));
  started = (bool_t *) calloc(Nchunk, sizeof(bool_t));

  for (nt = 0;  nt < Nchunk;  nt++) {
    chunk[nt].work  = work;
    chunk[nt].arg   = arg;
    chunk[nt].nt    = nt;
    chunk[nt].begin = nt * size;
    chunk[nt].end   = (nt + 1) * size < N ? (nt + 1) * size : N;
  }
  for (nt = 1;  nt < Nchunk;  nt++) {
    if (pthread_create(&id[nt], &input.thread_attr, runChunk, &chunk[nt])) {

      /* --- Could not start the thread: do its share here instead --- */

      Error(WARNING, "parallelFor",
	    "Could not create a thread, running its share in the caller");
      runChunk(&chunk[nt]);
    } else
      started[nt] = TRUE;
  }
  runChunk(&chunk[0]);
  for (nt = 1;  nt < Nchunk;  nt++)
    if (started[nt]) pthread_join(id[nt], NULL);

  free(started);
  free(id);
  free(chunk);
}
