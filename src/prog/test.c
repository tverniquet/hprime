#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>

#include "misc.h"
#include "prime_count.h"

int
main (int argc, char *argv[])
{
   uint64_t max;
   uint64_t s;
   uint64_t count;
   int plan_0;
   int plan_1;
   int nthreads = 0;

   if (argc < 5)
      exit_error("usage: %s min max plan_0 plan_1 nthreads\n", argv[0]);


   s = strtol(argv[1], NULL, 0);
   max = strtol(argv[2], NULL, 0);
   plan_0 = atoi(argv[3]);
   plan_1 = atoi(argv[4]);
   if (argc > 5)
      nthreads = atoi(argv[5]);

   getprimecount_cmp_plan(s, max, &count, plan_0, plan_1, nthreads);

   printf ("%"PRIu64"\n", count);

   return EXIT_SUCCESS;
}
