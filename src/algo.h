/*

	This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/

#ifndef _ALGO_H_
#define _ALGO_H_

#include<vector>
#include<map>

#include "globals.h"
#include "taxonomy.h"
#include "utility.h"

struct pathNode_st{
	IDnum taxonID;
	float likelihood;
	int category; // 0->root, 1->phylum, 2->genus, 3->species;
	pathNode_st *prevNode;
};

struct sequence_st{
	string seqName;
	vector<Gene> genes;
	map<IDnum, PathNode*> seqTaxonForest;
	
	void printSeq();
};

struct gene_st{
	vector<IDnum> gis;
	vector<float> identity;
	vector<float> bitscore;
	vector<IDnum> clusters;
	vector<float> dualHist;
	vector<float> subMTX;
	vector<IDnum> taxonIDs;
};

// initializers and destroyers;
Sequence *newSequence();
Gene *newGene();
PathNode *newPathNode();

vector<string> split(string s, char delim);

// load information from input file
vector<Sequence> loadInfoFromInputFile(const char* infile, int N);

void loadGI2TaxonLibFromFile(const char* gi2taxonFile, vector<Sequence> &QuerySeq);

void loadGI2ClstrLibFromFile(const char* gi2clstrFile, vector<Sequence> &QuerySeq);

void likelihoodCal(TaxonTree *tTree, vector<Sequence> &QuerySeq);

void writeResultsToOutputFile(const char* outfile, TaxonTree *tTree, TaxonName *tName,
								 vector<Sequence> &QuerySeq, float thr);

#endif