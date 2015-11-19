#ifndef _HARU_MISC_H
#define _HARU_MISC_H

#include <inttypes.h>


/* Obviously only useful by directly specifying a static array */
#define ARR_SIZEOF(A) (sizeof A / sizeof A[0])

/*
 * Macro allows use for constant assignment. These suffer side effects so use
 * carefully.
 */
#define CEIL_DIV(A, D) (((A) + (D)-1) / (D))
#define CEIL_TO(A, D)  (CEIL_DIV(A, D) * (D))

#define FLOOR_TO(A, D)  ((A) / (D) * (D))


#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))


#define FREE(A) do if (A) { free(A); (A) = NULL; } while (0)


/*
 * For vector stuff
 */
typedef unsigned char v16ui   __attribute__((vector_size(16)));
typedef uint16_t      v16_8ui __attribute__((vector_size(16)));


typedef unsigned char v32ui   __attribute__((vector_size(32)));
typedef uint64_t      v32_4ui __attribute__((vector_size(32)));
typedef int64_t       v32_4si __attribute__((vector_size(32)));
typedef uint32_t      v32_8ui __attribute__((vector_size(32)));
typedef int32_t       v32_8si __attribute__((vector_size(32)));
typedef uint16_t      v32_16ui __attribute__((vector_size(32)));
typedef int16_t       v32_16si __attribute__((vector_size(32)));

void print_v16 (v16ui *a);
void print_v32 (v32ui *a);
void print_v32_4ui (v32_4ui *a);



/****************************************************************************
 * Some extra funky GCC macros.
 *
 * Hint: To view the expansion compile with -E (and run indent)
 *
 * The important macro at the bottom is the "DO_FOR", it essentially expands:
 *
 * #define FOO(a) printf("%d\n", a);
 *
 * DO_FOR(FOO, 6,3,1,3,5,6)
 *
 * to:
 * printf("%d\n", 6);
 * printf("%d\n", 3);
 * printf("%d\n", 1);
 * printf("%d\n", 3);
 * printf("%d\n", 5);
 * printf("%d\n", 6);
 *
 *
 * The main point is to force compile time constants. I noticed the compiler
 * often had difficulty resolving/using the immediate value otherwise.
 *
 ****************************************************************************/
#define  _NARGS2(X,M,N18,N17,N16,N15,N14,N13,N12,N11,N10,N9,N8,N7,N6,N5,N4,N3,N2,N1,N, ...) N
#define NA(...) _NARGS2(0, 0, __VA_ARGS__, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#define  L1(N, M, VAL,...) M(VAL)
#define  L2(N, M, VAL,...) M(VAL) L1(N, M, __VA_ARGS__)
#define  L3(N, M, VAL,...) M(VAL) L2(N, M, __VA_ARGS__)
#define  L4(N, M, VAL,...) M(VAL) L3(N, M, __VA_ARGS__)
#define  L5(N, M, VAL,...) M(VAL) L4(N, M, __VA_ARGS__)
#define  L6(N, M, VAL,...) M(VAL) L5(N, M, __VA_ARGS__)
#define  L7(N, M, VAL,...) M(VAL) L6(N, M, __VA_ARGS__)
#define  L8(N, M, VAL,...) M(VAL) L7(N, M, __VA_ARGS__)
#define  L9(N, M, VAL,...) M(VAL) L8(N, M, __VA_ARGS__)
#define L10(N, M, VAL,...) M(VAL) L9(N, M, __VA_ARGS__)
#define L11(N, M, VAL,...) M(VAL) L10(N, M, __VA_ARGS__)
#define L12(N, M, VAL,...) M(VAL) L11(N, M, __VA_ARGS__)
#define L13(N, M, VAL,...) M(VAL) L12(N, M, __VA_ARGS__)
#define L14(N, M, VAL,...) M(VAL) L13(N, M, __VA_ARGS__)
#define L15(N, M, VAL,...) M(VAL) L14(N, M, __VA_ARGS__)
#define L16(N, M, VAL,...) M(VAL) L15(N, M, __VA_ARGS__)
#define L17(N, M, VAL,...) M(VAL) L16(N, M, __VA_ARGS__)
#define L18(N, M, VAL,...) M(VAL) L17(N, M, __VA_ARGS__)

#define LX(N, M, ...)  L##N(N, M, __VA_ARGS__)
#define LXN(N, M, ...) LX(N, M, __VA_ARGS__)

#define DO_FOR(M, ...) LXN(NA(__VA_ARGS__), M, __VA_ARGS__)
/****************************************************************************/


/*
 * NOTE:
 * I use clang as the syntax checker but compile with gcc. But clang doesn't
 * seem to know about shuffle. This just stops __builtin_shuffle showing an
 * error in my syntax checker.
 */
#ifdef __clang__
#define DO_NONE(A) (void)A,
#define __builtin_shuffle(A, ...) (DO_FOR(DO_NONE, __VA_ARGS__) A)
#define __builtin_assume_aligned(A,B) A
#endif


/*
 * Used for timing various things
 */
#include <time.h>
#include <inttypes.h>

#define TIME_DIFF_FMT_MS    "[%"PRIu64"ms, %"PRIu64"ms]"
#define TIME_DIFF_VALUES_MS(TS) (TS)->c_diff/1000, (TS)->r_diff/1000

#define TIME_DIFF_FMT_US    "[%"PRIu64"us, %"PRIu64"us]"
#define TIME_DIFF_VALUES_US(TS) (TS)->c_diff, (TS)->r_diff

struct timespot {
   clock_t c;
   struct timespec r;
   uint64_t c_diff; /* usecs */
   uint64_t r_diff; /* usecs */
};


/*
 * Init the timespot (bzero it)
 */
void init_time(struct timespot *t);


/*
 * Mark the start of a session
 */
void mark_time(struct timespot *t);


/*
 * Add the time difference from the last mark_time() call to the cumulated time
 * score for the timespot.
 */
void add_timediff(struct timespot *t);


/*
 * Init the timespot (bzero it)
 */
void get_time_string(struct timespot *t);


/*
 * Program related
 */

#define exit_error(...) \
   do { \
      fprintf(stderr, __VA_ARGS__); \
      fprintf(stderr, "%d:%s:%s\n", __LINE__, __FILE__, __func__); \
      exit(1); \
   } while (0)



#endif /* _HARU_MISC_H */
