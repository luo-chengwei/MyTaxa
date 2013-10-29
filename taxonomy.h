/*

	This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/

#ifndef _TAXONOMY_H_
#define _TAXONOMY_H_

#include <map>
#include <vector>
#include "globals.h"

using namespace std;

struct taxonNode_st {
	IDnum taxonID;
	taxonNode_st *prevNode;
	char *rank;
};

struct nameRank_st {
	string name;
	string rank;
};

struct IDRank_st{
	IDnum taxonID;
	string rank;
};

struct taxonTree_st {
	map<IDnum, TaxonNode*> nodes;
};

struct taxonName_st {
	map<IDnum, string> names;
};

// initializer and destroyer
TaxonNode *newTaxonNode();

TaxonTree *newTaxonTree();

void destroyTaxonTree(TaxonTree *tTree);

TaxonName *newTaxonName();

void destroyTaxonName(TaxonName *tName);

// load database from db files;
TaxonTree *importTaxonTreeFromFile(const char* taxonTreeFile);

TaxonName *importTaxonNameFromFile(const char* taxonNameFile);

// utility functions that are useful in runtime
vector<NameRank> taxonomyPath(TaxonTree *tTree, TaxonName *tNames, IDnum taxonID);

vector<IDRank> taxonomyPathIDRank(TaxonTree *tTree, IDnum taxonID);

vector<IDnum> taxonomyPath(TaxonTree *tTree, IDnum taxonID);

string taxonomyPathString(vector<NameRank> path);

IDnum lowestCommonAncestor(TaxonTree *tTree, IDnum taxonIDA, IDnum taxonIDB);

#endif
