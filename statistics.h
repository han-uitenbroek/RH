/* ------- file: -------------------------- statistics.h ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jan 24 15:57:13 2000 --

       --------------------------                      ----------RH-- */

#ifndef __STATISTICS_H__
#define __STATISTICS_H__

/* --- Keep track of program statistics like CPU time and memory usage
       (not yet implemented).
       --                                              -------------- */

#define N_TIME_LEVELS  5

enum CPUaction  {TIME_START, TIME_POLL, TIME_ADD};

typedef struct {
  bool_t  printCPU;
  long    memory;
  FILE   *fp_CPU;
} ProgramStats;

/* --- Associated function prototypes --               -------------- */

void getCPU(int level, enum CPUaction action, char *label);
void printTotalCPU(void);

#endif /* !__STATISTICS_H__ */

/* ------- end ---------------------------- statistics.h ------------ */
