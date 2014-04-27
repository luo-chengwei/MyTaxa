/*
Copyright 2009 John Marshall (jm18@sanger.ac.uk) 

    This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>
#include <sys/time.h>

#include "globals.h"
#include "utility.h"

static void allocExitError(const char *function, unsigned long long count,
                          unsigned long long size, const char *name)
{
       if (size == 1)
               exitErrorf(EXIT_FAILURE, true,
                          "Can't %s %llu %ss",
                          function, count, name);
       else
               exitErrorf(EXIT_FAILURE, true,
                          "Can't %s %llu %ss totalling %llu bytes",
                          function, count, name, count * size);
}

void *mallocOrExit3(size_t count, size_t size, const char *name)
{
       void *p = malloc(count * size);
       if (p == NULL && count > 0)
               allocExitError("malloc", count, size, name);

       return p;
}

void *callocOrExit3(size_t count, size_t size, const char *name)
{
       void *p = calloc(count, size);
       if (p == NULL && count > 0)
               allocExitError("calloc", count, size, name);

       return p;
}

void *reallocOrExit4(void *ptr, size_t count, size_t size, const char *name)
{
       void *p = realloc(ptr, count * size);
       if (p == NULL && count > 0)
               allocExitError("realloc", count, size, name);

       return p;
}


static const char *programName = NULL;

void setProgramName(const char *name)
{
       programName = name;
}

void exitErrorf(int exitStatus, bool showErrno, const char *format, ...)
{
       int savedErrno = errno;
       va_list args;
       va_start(args, format);
       if (programName)
            fprintf(stderr, "%s: ", programName);
       vfprintf(stderr, format, args);
       if (showErrno){
        	fprintf(stderr, ": %s", strerror(savedErrno));
        }
       fprintf(stderr, "\n");
       va_end(args);
}


StringBuffer *newStringBuffer(size_t size)
{
	StringBuffer *buffer;

	buffer = callocOrExit(1, StringBuffer);
	if (size > 0)
	{
		buffer->str = callocOrExit(size, char);
		buffer->allocated = size;
	}
	buffer->length = 1; /* virtual count for the final '\0' */

	return buffer;
}

void destroyStringBuffer(StringBuffer *buffer, bool freeString)
{
	if (!buffer)
		return;

	if (freeString && buffer->allocated > 0)
		free(buffer->str);
	free(buffer);
}

void appendStringBuffer(StringBuffer *buffer, char *str)
{
	int strSize;
	int newSize;

	if (!buffer)
		return;

	strSize = strlen(str);
	newSize = buffer->allocated;

	while (strSize + (int) buffer->length > newSize)
		newSize *= 2;
	if (newSize != (int) buffer->allocated)
	{
		buffer->str = reallocOrExit(buffer->str, newSize, char);
		if (!buffer->allocated)
			*buffer->str = '\0';
		buffer->allocated = newSize;
	}
	buffer->length += strSize;
	buffer->str = strcat(buffer->str, str);
}

void resetStringBuffer(StringBuffer *buffer)
{
	if (buffer && buffer->allocated > 0)
	{
		buffer->length = 1; /* Virtual count for the final '\0' */
		*buffer->str = '\0';
	}
}
