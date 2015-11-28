#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> /* XXX Debug */

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct prime_and_n
{
   uint16_t prime;
   uint16_t n;
}__attribute__((__packed__));


struct simple_middle_thread_ctx
{
   uint16_t *offsets;
   int       calculated_index;
   int64_t   last_blockno;
};


struct simple_middle_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   int      nthreads;
   int      primelist_count;
   struct prime_and_n              *primelist;
   struct simple_middle_thread_ctx *thread_data;
};


int
simple_middle_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   int i;
   struct simple_middle_ctx *sctx = malloc(sizeof (struct simple_middle_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime <= 1<<16);

   sctx->nthreads = pctx->num_threads ?: 1;
   sctx->block_size = pctx->current_block->block_size;
   sctx->thread_data = malloc(sctx->nthreads * sizeof(struct simple_middle_thread_ctx));

   sctx->primelist = malloc(sizeof(struct prime_and_n) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;

   for (i = 0; i < sctx->nthreads; i++) {
      sctx->thread_data[i].offsets = malloc(sizeof(uint16_t ) * 8 * (sctx->end_prime - start_prime) / 4 + 1000);
      sctx->thread_data[i].calculated_index = 0;
      sctx->thread_data[i].last_blockno = INT64_MAX;
   }
   return 0;
}


int
simple_middle_free(void *ctx)
{
   int i;
   struct simple_middle_ctx *sctx = ctx;

   for(i = 0; i < sctx->nthreads; i++)
      FREE(sctx->thread_data[i].offsets);

   FREE(sctx->primelist);
   FREE(sctx->thread_data);
   FREE(ctx);
   return 0;
}


int
simple_middle_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct simple_middle_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      sctx->primelist[sctx->primelist_count].prime = primelist[*ind];
      sctx->primelist[sctx->primelist_count].n = sctx->block_size / primelist[*ind];
      sctx->primelist_count++;
   }
   return 0;
}


int
simple_middle_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx)
{
   (void)pctx;
   (void)target_num;

   struct simple_middle_ctx *sctx = ctx;
   /* Recalculate all of the offsets depending on the new block */
   sctx->thread_data[pctx->thread_index].calculated_index = 0;
   return 0;
}


static void
compute_block_first_time(struct prime_current_block *pcb, uint32_t sieve_prime, uint16_t *next_offsets)
{
   char     *bmp;
   int       i;
   int       offsets[8];
   const unsigned char bits[]  = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};

   for (i = 0; i < 8; i++)
      offsets[pp_to_bit(ind_to_mod[i] * sieve_prime)] = num_to_bytes(ind_to_mod[i] * sieve_prime);

   bmp = pcb_initial_offset(pcb, sieve_prime);

   for (i = 0; i < 8; i++)
      if (pcb_inrange(pcb, bmp + offsets[i]))
         *(bmp + offsets[i]) |= bits[i];

   bmp += (bmp + sieve_prime < pcb_end(pcb)) ? sieve_prime : 0;

   for ( ; bmp < pcb_end(pcb) - sieve_prime; bmp += sieve_prime)
      for (i = 0; i < 8; i++)
         *(bmp + offsets[i]) |= bits[i];

   for (i = 0; i < 8; i++)
      if (pcb_inrange(pcb, bmp + offsets[i]))
         *(bmp + offsets[i]) |= bits[i];

   for (i = 0; i < 8; i++)
      next_offsets[i] = (bmp + offsets[i]) - pcb_end(pcb) + ((bmp + offsets[i] < pcb_end(pcb)) ? sieve_prime : 0);
}


static void
check_new_sieve_primes(struct simple_middle_ctx *sctx, struct simple_middle_thread_ctx *tdata, struct prime_current_block *pcb)
{
   uint32_t sieve_prime;

   for ( ; tdata->calculated_index < sctx->primelist_count; tdata->calculated_index++) {

      sieve_prime = sctx->primelist[tdata->calculated_index].prime;

      if (sieve_prime > pcb->sqrt_end_num)
         return;

      compute_block_first_time(pcb, sctx->primelist[tdata->calculated_index].prime, &tdata->offsets[tdata->calculated_index*8]);
   }
}


static void
compute_block(struct prime_current_block *pcb, uint32_t sieve_prime, uint16_t *offsets, uint32_t n)
{
   char *bmp;
   int i;
   const unsigned char bits[] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};

   for (bmp = pcb->block; n--; bmp += sieve_prime)
      for (i = 0; i < 8; i++)
         *(bmp + offsets[i]) |= bits[i];

   for (i = 0; i < 8; i++) {
      if ((bmp + offsets[i]) < pcb_end(pcb))
         *(bmp + offsets[i]) |= bits[i];
      offsets[i] = (bmp + offsets[i]) - pcb_end(pcb) + ((bmp + offsets[i]) < pcb_end(pcb) ? sieve_prime : 0);
   }
}


int
simple_middle_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct simple_middle_ctx        *sctx = ctx;
   struct simple_middle_thread_ctx *tdata = &sctx->thread_data[ptx->thread_index];
   /* Mark off multiples of 'a' in the block */
   int  i;
   int64_t  skip = ptx->current_block.block_num - tdata->last_blockno;

   tdata->last_blockno = ptx->current_block.block_num;

   if (sctx->primelist_count == 0)
      return 0;

   if (sctx->primelist[0].prime == 7)
      memset(ptx->current_block.block, 0, ptx->current_block.block_size);

   if (skip < 0 || skip > 1) {
      tdata->calculated_index = 0;
      skip = 0;
   }

   for (i = 0; i < tdata->calculated_index; i++)
      compute_block(&ptx->current_block, sctx->primelist[i].prime, &tdata->offsets[i*8], sctx->primelist[i].n);

   check_new_sieve_primes(sctx, tdata, &ptx->current_block);

   return 0;
}
