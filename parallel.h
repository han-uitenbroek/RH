/* ------- file: -------------------------- parallel.h --------------
       --------------------------                      ----------RH-- */

#ifndef __RH_PARALLEL_H__
#define __RH_PARALLEL_H__

/* --- Work on the half-open range [begin, end) of items. nt is the
       index of the worker, in 0..nthreads-1.  --      -------------- */

typedef void (*ParallelRange)(int begin, int end, int nt, void *arg);

void parallelFor(int N, int nthreads, ParallelRange work, void *arg);

#endif /* __RH_PARALLEL_H__ */
