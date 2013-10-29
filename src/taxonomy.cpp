/*

	This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/

#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <vector>

#include "taxonomy.h"
#include "utility.h"
#include "globals.h"
#include "algo.h"

using namespace std;


// initializers and destroyers

TaxonNode *newTaxonNode(){
	TaxonNode *tNode = callocOrExit(1, TaxonNode);
	return tNode;
}

TaxonTree *newTaxonTree(){
	TaxonTree *tTree = callocOrExit(1, TaxonTree);
	tTree->nodes.clear();
	return tTree;
}

void destroyTaxonTree(TaxonTree *tTree){
	if(tTree == NULL){
		return;
	}
	
	if(tTree->nodes.size() != 0){
		tTree->nodes.clear();
	}
	
	free(tTree);
}

TaxonName *newTaxonName(){
	TaxonName *tName = callocOrExit(1, TaxonName);
	tName->names.clear();
	return tName;
}

void destroyTaxonName(TaxonName *tName){
	if(tName == NULL){
		return;
	}
	
	if(tName->names.size() != 0){
		tName->names.clear();
	}
	
	free(tName);
}

TaxonNode *taxonNodeInTaxonTree(TaxonTree *tTree, IDnum taxonID){
	try{
		return tTree->nodes[taxonID];
	}catch(...){
		return NULL;
	}
}

void addNodeToTaxonTree(TaxonTree *tTree, IDnum nodeIDnum, 
				IDnum prevNodeIDnum, char *rank){
	
	TaxonNode *currentNode = taxonNodeInTaxonTree(tTree, nodeIDnum);
	TaxonNode *prevNode = taxonNodeInTaxonTree(tTree, prevNodeIDnum);
	
	if(currentNode == NULL){
		currentNode = newTaxonNode();
		currentNode->taxonID = nodeIDnum;
		currentNode->rank = rank;
		tTree->nodes[nodeIDnum]=currentNode;	
	}

	
	if(prevNode == NULL){
		prevNode = newTaxonNode();
		prevNode->taxonID = prevNodeIDnum;
		tTree->nodes[prevNodeIDnum]=prevNode;
	}
	
	if(nodeIDnum != prevNodeIDnum){
		tTree->nodes[nodeIDnum]->prevNode=tTree->nodes[prevNodeIDnum];
	}
	
	if(currentNode->rank == NULL){
		currentNode->rank = rank;
	}
}

// function that reads taxonNodes lib from NCBI file
TaxonTree *importTaxonTreeFromFile(const char* taxonTreeFile){
	FILE *ncbiTaxonTreeFile = fopen(taxonTreeFile, "r");
	
	const int maxLine = 20000;
	char line[maxLine];
	IDnum currentNode;
	IDnum prevNode;
	char *rank;
	
	TaxonTree *tTree = newTaxonTree();
	
//	cout << "Now reading NCBI taxonomy file: " << taxonFile << endl;
	if (ncbiTaxonTreeFile == NULL){
		cerr << "Could not open input NCBI taxonomy file: " << taxonTreeFile << endl;
		exit(EXIT_FAILURE);
	}
	
	while(fgets(line, maxLine, ncbiTaxonTreeFile) != NULL){
		string tmpA, tmpB, tmpRank;		
		stringstream streamLine;
		streamLine << line;
		try{
			streamLine >> currentNode >> prevNode >> tmpA >> tmpB;
			if(tmpB.empty()){
				tmpRank = tmpA;
			}else{
				tmpRank = tmpA + string(" ") + tmpB;
			}
			rank = new char[tmpRank.size()+1];
			copy(tmpRank.begin(), tmpRank.end(), rank);
			rank[tmpRank.size()] = '\0';
		}catch(...){
			streamLine >> currentNode >> prevNode >> tmpA;
			rank = new char[tmpA.size()];
			copy(tmpA.begin(), tmpA.end(), rank);
			rank[tmpA.size()] = '\0';
		}

		addNodeToTaxonTree(tTree, currentNode, prevNode, rank);
	}
	
	fclose(ncbiTaxonTreeFile);
	
	return tTree;
}

TaxonName *importTaxonNameFromFile(const char* taxonNameFile){
	FILE *ncbiTaxonNameFile = fopen(taxonNameFile, "r");
	
	const int maxLine = 5000;
	char line[maxLine];
	IDnum taxonID;
	string taxonName;
	char delim = '\t';
	TaxonName *tName = newTaxonName();
	
	if (ncbiTaxonNameFile == NULL){
		cerr << "Could not open input NCBI taxonomy file: " << taxonNameFile << endl;
		exit(EXIT_FAILURE);
	}
	
	
	while(fgets(line, maxLine, ncbiTaxonNameFile) != NULL){
		vector<string> elems = split(string(line), delim);
		taxonName = elems[1];
		taxonID = atoi(elems[0].c_str());
		
		map<IDnum, string>::iterator it;
		it = tName->names.begin();
		tName->names.insert(it, pair<IDnum, string>(taxonID, taxonName));
	}
	
	fclose(ncbiTaxonNameFile);
	
	return tName;
}

// some operational functions
bool isRoot(TaxonNode* tNode){
	if(tNode->taxonID == 1){
		return true;
	}else{
		return false;
	}
}

// output the taxonomy path in a vector<NameRank>, given a taxonID;
vector<NameRank> taxonomyPath(TaxonTree *tTree, TaxonName *tNames, IDnum taxonID){
	vector<NameRank> taxonPath;
	NameRank nr;
	TaxonNode *startNode = newTaxonNode();
	TaxonNode *currentNode = newTaxonNode();
	
	try{
		startNode = tTree->nodes[taxonID];
	}catch(...){
		cout << "taxonID: " << taxonID << " not found in database"<< endl;
	}
	
	// visit every node from leaf to root;
	try{
		currentNode = startNode; 
		while(!isRoot(currentNode)){
			IDnum currentTaxonID = currentNode->taxonID;
			string name = tNames->names[currentTaxonID];
			char *rank = tTree->nodes[currentTaxonID]->rank;
			nr.name = name;
			nr.rank = string(rank);
			taxonPath.push_back(nr);
			currentNode = currentNode->prevNode;
		}
	}catch(...){
		cout << "taxonID : " << currentNode->taxonID << " not found"<<endl;
	}
	
	return taxonPath;
}

// reload returns vector<IDRank>;
vector<IDRank> taxonomyPathIDRank(TaxonTree *tTree, IDnum taxonID){
	vector<IDRank> taxonPath;
	IDRank idr;
	TaxonNode *currentNode = newTaxonNode();
	
	try{
		currentNode = tTree->nodes[taxonID];
	}catch(...){
		cout << "taxonID: " << taxonID << " not found in database"<< endl;
	}
	
	
	// visit every node from leaf to root;
	if(currentNode == NULL){
		return taxonPath;
	}
	
	try{
		while(!isRoot(currentNode)){
			IDnum currentTaxonID = currentNode->taxonID;
			char *rank = tTree->nodes[currentTaxonID]->rank;
			idr.taxonID = currentTaxonID;
			idr.rank = string(rank);
			taxonPath.push_back(idr);
			currentNode = currentNode->prevNode;
		}
	}catch(...){
		cout << "taxonID : " << currentNode->taxonID << " not found"<<endl;
	}
	
	return taxonPath;
}

// reload taxonomyPath, returns vector<IDnum>;
vector<IDnum> taxonomyPath(TaxonTree *tTree, IDnum taxonID){
	vector<IDnum> taxonPath;
	
	TaxonNode *startNode = newTaxonNode();
	TaxonNode *currentNode = newTaxonNode();
	
	try{
		startNode = tTree->nodes[taxonID];
	}catch(...){
		cout << "taxonID: " << taxonID << " not found in database"<< endl;
	}
	
	// visit every node from leaf to root;
	try{
		currentNode = startNode; 
		while(!isRoot(currentNode)){
			IDnum currentTaxonID = currentNode->taxonID;
			taxonPath.push_back(currentTaxonID);
			currentNode = currentNode->prevNode;
		}
	}catch(...){
		cout << "taxonID : " << currentNode->taxonID << " not found"<<endl;
	}
	
	return taxonPath;
}

string taxonomyPathString(vector<NameRank> path){
	string pathString, tmpString;
	
	for(int index = path.size()-1; index >= 0; index--){
		NameRank node = path[index];
		if(node.rank.find("no rank") == string::npos){
			tmpString += "<" + node.rank + ">" + node.name.substr(0, node.name.size()-1) + ";";
		}else if(node.rank.find("group") != string::npos){
			tmpString += node.name.substr(0, node.name.size()-1) + ";";
		}
	}
	
	pathString = tmpString.substr(0, tmpString.size()-1);
	
	return pathString;
}


// return the lowest common ancestor (LCA) of two taxonIDs.
IDnum lowestCommonAncestor(TaxonTree *tTree, IDnum taxonIDA, IDnum taxonIDB){
	vector<IDnum> taxonPathA = taxonomyPath(tTree, taxonIDA);
	vector<IDnum> taxonPathB = taxonomyPath(tTree, taxonIDB);
	
	vector<IDnum>::iterator iterA;
	vector<IDnum>::iterator iterB;
	
	for(iterA = taxonPathA.begin(); iterA != taxonPathA.end(); iterA++){
		for(iterB = taxonPathB.begin(); iterB != taxonPathB.end(); iterB++){
			if(*iterA == *iterB){
				return *iterA;
			}
		}
	}
	
	return 1;
}



