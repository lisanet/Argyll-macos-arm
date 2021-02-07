/*
 * Copyright (c) 2007-2014, Lloyd Hilaiel <me@lloyd.io>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#ifndef __YAJL_COMMON_H__
#define __YAJL_COMMON_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define YAJL_MAX_DEPTH 128

// We're not creating a DLL, so don't mark the API's - GWG

#ifdef NEVER

/* msft dll export gunk.  To build a DLL on windows, you
 * must define WIN32, YAJL_SHARED, and YAJL_BUILD.  To use a shared
 * DLL, you must define YAJL_SHARED and WIN32 */
#if (defined(_WIN32) || defined(WIN32)) && defined(YAJL_SHARED)
#  ifdef YAJL_BUILD
#    define YAJL_API __declspec(dllexport)
#  else
#    define YAJL_API __declspec(dllimport)
#  endif
#else
#  if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 303
#    define YAJL_API __attribute__ ((visibility("default")))
#  else
#    define YAJL_API
#  endif
#endif

#else

# define YAJL_API

#endif

// Create a cross platform 64 bit int type "longlong" - GWG
#ifndef NUMLIB_H

#if (__STDC_VERSION__ >= 199901L)	/* C99 */

#include <stdint.h> 

typedef int64_t longlong ;

#define PF64PREC "ll"		/* printf format precision specifier */
#define CF64PREC "LL"		/* Constant precision specifier */

#ifndef LLONG_MIN
# define LLONG_MIN INT64_MIN
#endif

#ifndef LLONG_MAX
# define LLONG_MAX INT64_MAX
#endif

#else  /* !__STDC_VERSION__ */
#ifdef _MSC_VER

typedef __int64 longlong;

#define PF64PREC "I64"				/* printf format precision specifier */
#define CF64PREC "LL"				/* Constant precision specifier */

#ifndef LLONG_MIN
# define LLONG_MIN _I64_MIN
#endif

#ifndef LLONG_MAX
# define LLONG_MAX _UI64_MAX
#endif

#else  /* !_MSC_VER */

/* The following works on a lot of modern systems, including */
/* LLP64 and LP64 models, but won't work with ILP64 which needs int32 */

#ifdef __GNUC__
# ifdef __LP64__	/* long long could be 128 bit ? */
   typedef long longlong;
#  define PF64PREC "l"			/* printf format precision specifier */
#  define CF64PREC "L"			/* Constant precision specifier */
#  ifndef LLONG_MAX
#   define LLONG_MAX	__LONG_MAX__
#  endif
#  ifndef LLONG_MIN
#   define LLONG_MIN	(-LLONG_MAX-1)
#  endif
#  ifndef ULLONG_MAX
#   define ULLONG_MAX	(LLONG_MAX * 2UL + 1)
#  endif
# else				/* long could be 32 bits */
   typedef long long longlong;
#  define PF64PREC "ll"			/* printf format precision specifier */
#  define CF64PREC "LL"			/* Constant precision specifier */
#  ifndef LLONG_MAX
#   define LLONG_MAX	__LONG_LONG_MAX__
#  endif
#  ifndef LLONG_MIN
#   define LLONG_MIN	(-LLONG_MAX-1)
#  endif
#  ifndef ULLONG_MAX
#   define ULLONG_MAX	(LLONG_MAX * 2ULL + 1)
#  endif
# endif /* !__LP64__ */
#endif /* __GNUC__ */

#endif	/* !_MSC_VER */
#endif  /* !__STDC_VERSION__ */

#else /* !NUMLIB_H */

typedef INR64 longlong ;

#endif /* !NUMLIB_H */

/** pointer to a malloc function, supporting client overriding memory
 *  allocation routines */
typedef void * (*yajl_malloc_func)(void *ctx, size_t sz);

/** pointer to a free function, supporting client overriding memory
 *  allocation routines */
typedef void (*yajl_free_func)(void *ctx, void * ptr);

/** pointer to a realloc function which can resize an allocation. */
typedef void * (*yajl_realloc_func)(void *ctx, void * ptr, size_t sz);

/** A structure which can be passed to yajl_*_alloc routines to allow the
 *  client to specify memory allocation functions to be used. */
typedef struct
{
    /** pointer to a function that can allocate uninitialized memory */
    yajl_malloc_func malloc;
    /** pointer to a function that can resize memory allocations */
    yajl_realloc_func realloc;
    /** pointer to a function that can free memory allocated using
     *  reallocFunction or mallocFunction */
    yajl_free_func free;
    /** a context pointer that will be passed to above allocation routines */
    void * ctx;
} yajl_alloc_funcs;

#ifdef __cplusplus
}
#endif

#endif
