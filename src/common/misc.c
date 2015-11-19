#include <string.h>
#include <stdio.h>


#include "misc.h"

void
init_time (struct timespot *t)
{
   bzero(t, sizeof *t);
}


void
mark_time (struct timespot *t)
{
   t->c = clock();
   clock_gettime(CLOCK_MONOTONIC, &t->r);
}


void
add_timediff (struct timespot *t)
{
   struct timespot l;
   mark_time(&l);
   t->c_diff += (l.c - t->c)*1000000 / CLOCKS_PER_SEC;
   t->r_diff += (l.r.tv_sec - t->r.tv_sec) * 1000000 + ((long)l.r.tv_nsec - t->r.tv_nsec) / 1000;
}


void
print_v16 (v16ui *a)
{
   unsigned char t2[16] __attribute__((aligned(16)));
   *(v16ui *)t2 = *a;
   int i;
   for (i = 0; i < 16; i++) {
      fprintf(stderr, "%02X ", t2[i]);
      if (i%4 == 3)
         fputc(' ', stderr);
   }
   fputs("\n", stderr);
}


void
print_v32(v32ui *a)
{
   unsigned char t2[32] __attribute__((aligned(32)));
   *(v32ui *)t2 = *a;
   int i;
   for (i = 0; i < 32; i++) {
      fprintf(stderr, "%02X ", t2[i]);
      if (i%4 == 3)
         fputc(' ', stderr);
   }
   fputs("\n", stderr);
}


void
print_v32_4ui(v32_4ui *a)
{
   uint64_t b[4] __attribute__((aligned(32)));
   int i;
   *(v32_4ui *)b = *a;
   for (i = 0; i < 4; i++) {
      fprintf(stderr, "%6lu ", b[i]);
   }
}

