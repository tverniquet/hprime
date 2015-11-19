#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>

#include "prime_count.h"

#include "misc.h"

int
main (int argc, char *argv[])
{
   uint64_t max;
   uint64_t s;
   uint64_t count;
   int ind = 0;
   int nthreads = 0;
   int inorder = 0;
   struct timespot ts;

   bzero(&ts, sizeof ts);

   if (argc < 3)
      exit_error("Usage: %s min max [plan] [nthreads] [inorder]\n", argv[0]);

   s = strtol(argv[1], NULL, 0);
   max = strtol(argv[2], NULL, 0);
   if (argc > 3)
      ind = strtol(argv[3], NULL, 0);

   if (argc > 4)
      nthreads = strtol(argv[4], NULL, 0);

   if (argc > 5)
      inorder = strtol(argv[5], NULL, 0);

   mark_time(&ts);
   getprimecount(ind, s, max, &count, nthreads, inorder);
   add_timediff(&ts);

   printf("%"PRIu64" "TIME_DIFF_FMT_MS"\n", count, TIME_DIFF_VALUES_MS(&ts));

   return EXIT_SUCCESS;
}
