#ifndef UTILITY_H
#define UTILITY_H

#include "TSP.h";

#define MALLOC(n, type) ((type *) malloc (n * sizeof(type)))
#define CALLOC(n, type) ((type *) calloc (n * sizeof(type)))
#define REALLOC(ptr, n, type) (realloc(ptr, n * sizeof(type)))

#endif
