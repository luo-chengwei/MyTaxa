/*
	Copyright 2009 John Marshall (jm18@sanger.ac.uk) 

    This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/
#ifndef UTILITY_H_
#define UTILITY_H_

#ifdef __GNUC__
#define ATTRIBUTE(list)  __attribute__ (list)
#else
#define ATTRIBUTE(list)
#endif

#include <stdio.h>

// Wrappers for malloc(), calloc(), and realloc() that always succeed.
// These functions print an error message and exit on failure, rather than
// requiring the calling function to check for NULL.  The arguments contain
// the type itself -- mallocOrExit(n, Foo) rather than malloc(n * sizeof Foo)
// -- to enable type checking and so that it can be shown in error messages.
#define mallocOrExit(count, type) \
               ((type *) mallocOrExit3((count), sizeof(type), #type))
#define callocOrExit(count, type) \
               ((type *) callocOrExit3((count), sizeof(type), #type))
#define reallocOrExit(ptr, count, type) \
               ((type *) reallocOrExit4((ptr), (count), sizeof(type), #type))

// However there are types for which just appending a '*' produces the
// wrong type or a syntax error, rather than a pointer-to-<type>.  These
// less type-safe wrappers are provided for use in these unusual cases.
#define mallocOrExitWithoutCast(count, type) \
               (mallocOrExit3((count), sizeof(type), #type))
#define callocOrExitWithoutCast(count, type) \
               (callocOrExit3((count), sizeof(type), #type))
#define reallocOrExitWithoutCast(ptr, count, type) \
               (reallocOrExit4((ptr), (count), sizeof(type), #type))

// (Implementation functions -- use the macro wrappers above.)
void *mallocOrExit3(size_t count, size_t size, const char *name);
void *callocOrExit3(size_t count, size_t size, const char *name);
void *reallocOrExit4(void *ptr, size_t count, size_t size, const char *name);

// Sets the program name to be prepended to error messages.
void setProgramName(const char *name);

// Prints an error message to standard error (with printf-style formatting
// and optionally appending a perror-style description of errno), and calls
// exit() with the specified exit status.
void exitErrorf(int exitStatus, bool showErrno, const char *format, ...)
       ATTRIBUTE((format(printf, 3, 4), noreturn));


// String Buffer
typedef struct
{
	char *str;
	size_t length;
	size_t allocated;
}
StringBuffer;

StringBuffer *newStringBuffer(size_t size);
void destroyStringBuffer(StringBuffer *buffer, bool freeString);
void appendStringBuffer(StringBuffer *buffer, char *str);
void resetStringBuffer(StringBuffer *buffer);

#endif
