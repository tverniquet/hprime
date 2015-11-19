#include "wheel.h"

uint64_t
num_to_bytes(uint64_t num)
{
   return num / 30;
}


uint64_t
pp_to_bit(uint64_t num)
{

   /* Note: This only works by chance with an 8/30 wheel, otherwise a lookup
    * table is needed.
    *
    * Either might be faster.. there doesn't seem to be much difference between
    * the divide and using the lookup table (as in mod_to_ind).
    */

   return (num % 30)*8/30;
}


uint64_t
num_to_bit(uint64_t num)
{
   return mod_to_ind[num % 30];
}


uint64_t
bytes_to_num(uint64_t bytes)
{
   return bytes * 30;
}


uint64_t
get_next_prime(struct prime_current_block *pcb, uint32_t *byte, uint32_t *bit)
{
   uint64_t possible_prime_num = pcb->block_start_num + *byte * 30ul;
   unsigned char *bytep = (unsigned char *)pcb->block + *byte;

   for (;*byte < pcb->block_size;) {

      for (; *bit < 8; (*bit)++)
         if ((*bytep & (1ul<<*bit)) == 0)
            return possible_prime_num + (int)(ind_to_mod[(*bit)++] & 0xFF);

      possible_prime_num += 30;
      bytep++;
      (*byte)++;
      (*bit) = 0;
   }
   return 0;
}


uint64_t
get_next_possible_prime(uint64_t *byte, uint32_t *bit)
{
   if (++(*bit) == 8) {
      (*byte)++;
      *bit = 0;
   }
   return (*byte)*30 + ind_to_mod[*bit];
}


uint64_t
get_previous_possible_prime(uint64_t *byte, uint32_t *bit)
{
   if ((*bit)-- == 0) {
      (*byte)--;
      *bit = 7;
   }
   return (*byte)*30 + ind_to_mod[*bit];
}


const char skips[] = {1,6,4,2,4,2,4,6};

/*
 * The idea is to avoid needing to know 'b' and just incrementing the result
 *
 * a * b =   a * (b_byte)    ** input "base_byte"
 *         + a_bit * b_bit   ** calculated from tables **
 *         + a_byte * b_bit  ** increment b_bit
 *
 * The resulting base_byte will need to be ajusted via the a_x_b_bytes table
 */
void
get_next_result_set(uint32_t a_byte, uint32_t a_bit, uint64_t *base_byte, uint32_t *b_bit)
{
   if ((*b_bit) == 7) {
      *base_byte += a_byte + ind_to_mod[a_bit];
      *b_bit = 0;
   }
   else
      (*b_bit)++;

   switch(*b_bit) {
      case 1:
      case 7:
         *base_byte += a_byte;
         *base_byte += a_byte;
      case 2:
      case 4:
      case 6:
         *base_byte += a_byte;
         *base_byte += a_byte;
      case 3:
      case 5:
         *base_byte += a_byte;
      case 0:
         *base_byte += a_byte;
   };
}


/******************/




/* 4d]] :r!./src/gen_tables.pl */

const unsigned char a_x_b_bytes[8][8] = {
   {  0, 0, 0, 0, 0, 0, 0, 0 },
   {  0, 1, 2, 3, 3, 4, 5, 6 },
   {  0, 2, 4, 4, 6, 6, 8,10 },
   {  0, 3, 4, 5, 7, 8, 9,12 },
   {  0, 3, 6, 7, 9,10,13,16 },
   {  0, 4, 6, 8,10,12,14,18 },
   {  0, 5, 8, 9,13,14,17,22 },
   {  0, 6,10,12,16,18,22,28 }
};


const unsigned char a_x_b_byte_diffs[8][8] = {
   {  0, 0, 0, 0, 0, 0, 0, 0 },
   {  0, 1, 1, 1, 0, 1, 1, 1 },
   {  0, 2, 2, 0, 2, 0, 2, 2 },
   {  0, 3, 1, 1, 2, 1, 1, 3 },
   {  0, 3, 3, 1, 2, 1, 3, 3 },
   {  0, 4, 2, 2, 2, 2, 2, 4 },
   {  0, 5, 3, 1, 4, 1, 3, 5 },
   {  0, 6, 4, 2, 4, 2, 4, 6 }
};


const unsigned char a_x_b_bits[8][8] = {
   {  0, 1, 2, 3, 4, 5, 6, 7 },
   {  1, 5, 4, 0, 7, 3, 2, 6 },
   {  2, 4, 0, 6, 1, 7, 3, 5 },
   {  3, 0, 6, 5, 2, 1, 7, 4 },
   {  4, 7, 1, 2, 5, 6, 0, 3 },
   {  5, 3, 7, 1, 6, 0, 4, 2 },
   {  6, 2, 3, 7, 0, 4, 5, 1 },
   {  7, 6, 5, 4, 3, 2, 1, 0 }
};


const unsigned char a_x_b_bitmask[8][8] = {
   {  1, 2, 4, 8,16,32,64,128 },
   {  2,32,16, 1,128, 8, 4,64 },
   {  4,16, 1,64, 2,128, 8,32 },
   {  8, 1,64,32, 4, 2,128,16 },
   { 16,128, 2, 4,32,64, 1, 8 },
   { 32, 8,128, 2,64, 1,16, 4 },
   { 64, 4, 8,128, 1,16,32, 2 },
   { 128,64,32,16, 8, 4, 2, 1 }
};


const unsigned char a_x_b_bitvals[8][8] = {
   {  1, 7,11,13,17,19,23,29 },
   {  7,19,17, 1,29,13,11,23 },
   { 11,17, 1,23, 7,29,13,19 },
   { 13, 1,23,19,11, 7,29,17 },
   { 17,29, 7,11,19,23, 1,13 },
   { 19,13,29, 7,23, 1,17,11 },
   { 23,11,13,29, 1,17,19, 7 },
   { 29,23,19,17,13,11, 7, 1 }
};


const unsigned int pp_diffs[8] = {6,4,2,4,2,4,6,2};

const char ind_to_mod[] = {1,7,11,13,17,19,23,29};

/* Note: rounds "up" to the next possible prime if not a "possible prime" */
const char mod_to_ind[] = {0,0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7};
