#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> /* XXX Debug */

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct prime_and_offset
{
   uint16_t prime;
   uint16_t offset;
}__attribute__((__packed__));


struct read_offs_thread_ctx {
   struct prime_and_offset *primelist;
   uint16_t *offsets;
   uint32_t primelist_count;
   uint32_t calculated_index;
   uint32_t top10_index[10];
   uint32_t top10_count[10];
   int      last_blockno;
};

struct read_offs_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   int nthreads;
   struct read_offs_thread_ctx *thread_data;
};


int
read_offs_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   int i;
   struct read_offs_ctx *sctx = calloc(1, sizeof (struct read_offs_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   sctx->nthreads = pctx->num_threads ?: 1;

   assert(sctx->end_prime <= pctx->current_block->block_size);
   assert(sctx->end_prime <= UINT16_MAX);

   sctx->thread_data = calloc(sizeof(struct read_offs_thread_ctx), sctx->nthreads);

   for (i = 0; i < sctx->nthreads; i++) {
      sctx->thread_data[i].primelist = malloc(sizeof(struct prime_and_offset) * (sctx->end_prime - start_prime) / 4 + 1000);
      sctx->thread_data[i].offsets = malloc(sizeof(uint16_t ) * 8 * (sctx->end_prime - start_prime) / 4 + 1000);
      sctx->thread_data[i].last_blockno = INT32_MAX;
      sctx->thread_data[i].calculated_index = 0;
   }
   return 0;
}


int
read_offs_free(void *ctx)
{
   int i;
   struct read_offs_ctx *sctx = ctx;
   for (i = 0; i < sctx->nthreads; i++) {
      FREE(sctx->thread_data[i].primelist);
      FREE(sctx->thread_data[i].offsets);
   }
   FREE(sctx->thread_data);
   FREE(ctx);
   return 0;
}


int
read_offs_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct read_offs_ctx *sctx = ctx;
   struct read_offs_thread_ctx *tdata;
   int i;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      for (i = 0; i < sctx->nthreads; i++) {
         tdata = &sctx->thread_data[i];
         tdata->primelist[tdata->primelist_count].prime = primelist[*ind];
         tdata->primelist[tdata->primelist_count].offset = primelist[*ind];
         tdata->primelist_count++;
      }
   }
   return 0;
}


int
read_offs_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx)
{
   (void)pctx;
   (void)target_num;

   struct read_offs_ctx *sctx = ctx;

   /* Recalculate all of the offsets depending on the new block */
   sctx->thread_data[pctx->thread_index].calculated_index = 0;
   sctx->thread_data[pctx->thread_index].last_blockno = INT32_MAX;

   return 0;
}


/*
 * The first time we don't know how many times to mark off the prime in the
 * block (ie don't really know 'n').
 */
static void
compute_block_first_time(struct prime_current_block *pcb, uint32_t sieve_prime, uint16_t *offsets)
{
   char *bmp;
   int   i;
   const unsigned char bits[]  = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};

   bmp = pcb->block + (uint64_t)sieve_prime * (sieve_prime / 30) - pcb->block_start_byte;

   for (; bmp < pcb->block + pcb->block_size; bmp += sieve_prime)
      for (i = 0; i < 8; i++)
         *(bmp + offsets[i]) |= bits[i];
}


static void
check_new_sieve_primes(struct read_offs_thread_ctx *tdata, struct prime_current_block *pcb, int skip, int mode)
{
   struct prime_and_offset *po;
   uint32_t sieve_prime;
   uint16_t *offsets;

   int i;

   for ( ; tdata->calculated_index < tdata->primelist_count; tdata->calculated_index++) {

      po = &tdata->primelist[tdata->calculated_index];
      sieve_prime = po->prime;
      offsets = &tdata->offsets[tdata->calculated_index * 8];

      for (i = 0; i < 8; i++)
         offsets[num_to_bit(ind_to_mod[i] * sieve_prime)] = num_to_bytes(ind_to_mod[i] * sieve_prime);

      /*
       * Mode 0 allows skipped primes to be added or reset primes to be reset and
       * still use the main mechanism (without using compute_block_first_time()
       */
      if (mode == 0) {
         if (sieve_prime * sieve_prime >= pcb->block_start_num)
            break;

         po->offset = sieve_prime - (pcb->block_start_byte - skip * pcb->block_size) % sieve_prime;
      }
      else {
         if (sieve_prime > pcb->sqrt_end_num)
            return;

         compute_block_first_time(pcb, sieve_prime, &tdata->offsets[tdata->calculated_index * 8]);

         po->offset = sieve_prime - pcb->block_start_byte % sieve_prime;
      }

      for (i = 10; i > 0; i--)
         if (sieve_prime <= pcb->block_size / i)
            tdata->top10_index[i-1] = tdata->calculated_index + 1;

      tdata->top10_count[9] = tdata->top10_index[9];
      for (i = 8; i >= 0; i--)
         if (tdata->top10_index[i] > 0)
            tdata->top10_count[i] = tdata->top10_index[i] - tdata->top10_index[i+1];

   }
}


static void
compute_get_set_offsets(struct prime_current_block *pcb, struct prime_and_offset *po, uint16_t *offsets, int *offs, int skip, int n)
{
   int32_t  sieve_prime = po->prime;
   int32_t  offset      = po->offset;
   int      i;

   while (skip--) {
      offset -= (pcb->block_size - n * sieve_prime);
      offset += (offset>>31) & sieve_prime;
   }

   po->offset = offset;
   offset -= sieve_prime;

   for (i = 0; i < 8; i++)
      offs[i] = (offset + offsets[i]) + ((int32_t)(offset + offsets[i])>>31 & sieve_prime);
}


/* the idea is that the inline remove the switch */
static inline void __attribute__((always_inline))
compute_block_n(struct prime_current_block *pcb, struct prime_and_offset *po, uint16_t *offsets, int skip, int n)
{
   const unsigned char bits[] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   char     *bmp;
   uint32_t  sieve_prime = po->prime;
   int       offs[8];
   int       i;

   compute_get_set_offsets(pcb, po, offsets, offs, skip, n - 1);

   bmp = pcb->block;

   switch (n) {
      case 10:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 9:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 8:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 7:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 6:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 5:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 4:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 3:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;
      case 2:
         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
         bmp += sieve_prime;

         for (i = 0; i < 8; i++)
            *(bmp + offs[i]) |= bits[i];
   }
}


static void
compute_block_low(struct prime_current_block *pcb, struct prime_and_offset *po, uint16_t *offsets, int skip)
{
   const unsigned char bits[] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   char    *bmp;
   int32_t  sieve_prime = po->prime;
   int      offs[8];
   int      n;
   int      i;

   n = pcb->block_size / sieve_prime + 1;

   compute_get_set_offsets(pcb, po, offsets, offs, skip, n - 1);

   for (bmp = pcb->block; n--; bmp += sieve_prime)
      for (i = 0; i < 8; i++)
         *(bmp + offs[i]) |= bits[i];
}



int
read_offs_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct read_offs_ctx        *sctx = ctx;
   struct read_offs_thread_ctx *tdata = &sctx->thread_data[ptx->thread_index];
   struct prime_and_offset     *po;
   int       skip = ptx->current_block.block_num - tdata->last_blockno;
   int       i;
   uint16_t *offs;

   if (tdata->primelist_count == 0)
      return 0;

   if (tdata->primelist[0].prime == 7)
      memset(ptx->current_block.block, 0, ptx->current_block.block_size);

   if (skip < 0 || skip > 8) {
      tdata->calculated_index = 0;
      skip = 0;
   }

   tdata->last_blockno = ptx->current_block.block_num;

   i = 0;

   offs = &tdata->offsets[0] - 8;
   po = tdata->primelist;

   check_new_sieve_primes(tdata, &ptx->current_block, skip, 0);

   /*
    * These didn't like being putting in a loop, presumably because of all of
    * the inlining, so I unrolled them manually
    */

   i = tdata->top10_count[9];
   while (i--)
      compute_block_low(&ptx->current_block, po++, (offs += 8), skip);

   i = tdata->top10_count[8];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 10);

   i = tdata->top10_count[7];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 9);

   i = tdata->top10_count[6];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 8);

   i = tdata->top10_count[5];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 7);

   i = tdata->top10_count[4];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 6);

   i = tdata->top10_count[3];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 5);

   i = tdata->top10_count[2];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 4);

   i = tdata->top10_count[1];
   while(i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 3);

   i = tdata->top10_count[0];
   while(i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), skip, 2);

   check_new_sieve_primes(tdata, &ptx->current_block, skip, 1);

   return 0;
}
