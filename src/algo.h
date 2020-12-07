/*

	This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/

#ifndef _ALGO_H_
#define _ALGO_H_

#include<vector>
#include<map>
#include<algorithm>

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
	
	float min_current_bitscore(){
		if(bitscore.size() == 0){
			return 0;
		}
		float min = 100000;
		for(vector<float>::iterator it = bitscore.begin(); it != bitscore.end(); ++it){
			if(min > *it){
				min = *it;
			}
		}
		return min;
	}
	
	float max_gap(){
		float max_gap = 0;
		for(unsigned int i = 0; i < bitscore.size() - 1; ++i){
			float min, max;
			min = (bitscore[i] < bitscore[i+1])?bitscore[i]:bitscore[i+1];
			max = (bitscore[i] > bitscore[i+1])?bitscore[i]:bitscore[i+1];
			float g = (max-min)/max;
			if(g > max_gap){
				max_gap = g;
			}
		}
		return max_gap;
	}
	
	void remove_min(){
		unsigned int min_index = 0;
		float min = 10000;
		for(vector<float>::iterator it = bitscore.begin(); it != bitscore.end(); ++it){
			if(min > *it){
				min = *it;
				min_index = it - bitscore.begin();
			}
		}
		gis.erase(min_index + gis.begin());
		bitscore.erase(min_index + bitscore.begin());
		identity.erase(min_index + identity.begin());
	}
	
	void organize_entries(){
		if(gis.size() < 2){
			return;
		}
		if(max_gap() > SCORE_DROP_THR){
			remove_min();
		}
	}
};

// initializers and destroyers;
Sequence *newSequence();
Gene *newGene();
PathNode *newPathNode();

vector<string> split(string s, char delim);

// load information from input file
vector<Sequence> loadInfoFromInputFile(const char* infile);

void loadGI2TaxonLibFromFile(const char* gi2taxonFile, vector<Sequence> &QuerySeq);

void loadGI2ClstrLibFromFile(const char* gi2clstrFile, vector<Sequence> &QuerySeq);

void likelihoodCal(TaxonTree *tTree, vector<Sequence> &QuerySeq);

void writeResultsToOutputFile(const char* outfile, TaxonTree *tTree, TaxonName *tName,
								 vector<Sequence> &QuerySeq, float thr);

#endif