/*

	This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/

#ifndef GLOBALS_H_
#define GLOBALS_H_

#ifndef NULL
#define NULL 0
#endif

#ifndef true
#define true 1
#define false 0
#endif

#ifndef __GNUC__
#define ATTRIBUTE_PACKED __attribute__ ((packed))
#else
#define ATTRIBUTE_PACKED
#endif

#define VERSION_NUMBER 1
#define RELEASE_NUMBER 0
#define UPDATE_NUMBER 0

// weights
#define W10  1
#define W11  1
#define W12  1
#define W20  1
#define W21  1
#define W22  1

// external structures here
struct taxonNode_st;
struct taxonTree_st;
struct taxonName_st;
struct nameRank_st;

// Namespace sizes here
#include <stdint.h>
typedef int32_t IDnum;
#endif

// Taxonomy elements
typedef struct taxonNode_st TaxonNode;
typedef struct taxonTree_st TaxonTree;
typedef struct taxonName_st TaxonName;
typedef struct nameRank_st NameRank;
typedef struct IDRank_st IDRank;

// algo elements
typedef struct sequence_st Sequence;
typedef struct gene_st Gene;
typedef struct pathNode_st PathNode;