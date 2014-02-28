/*

	This file is part of MeTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/


#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include "algo.h"
#include "utility.h"
#include "taxonomy.h"
#include "globals.h"

using namespace std;


Sequence *newSequence(){
	Sequence *seq = callocOrExit(1, Sequence);
	return seq;
}

Gene *newGene(){
	Gene *gene = callocOrExit(1, Gene);
	gene->gis.clear();
	gene->identity.clear();
	gene->bitscore.clear();
	gene->clusters.clear();
	gene->taxonIDs.clear();
	gene->dualHist.clear();
	gene->subMTX.clear();
	
	return gene;
}

PathNode *newPathNode(){
	PathNode *pNode = callocOrExit(1, PathNode);
	return pNode;
}

// get the Nth column of a given string (tab delimited);
vector<string> &split(string s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(string s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

// structure functions
void Sequence::printSeq(){
	cout << "Query sequence ID: " << seqName << endl;
	for(unsigned int i = 0; i < genes.size(); i++){
		cout << "  Gene " << i << endl;
		for(unsigned int j = 0; j < genes[i].gis.size(); j++){
			cout << " " << genes[i].gis[j] << " ";
			cout << " " << genes[i].identity[j] << " ";
			cout << " " << genes[i].bitscore[j] << " ";
			cout << " " << genes[i].taxonIDs[j] << " ";
			cout << " " << genes[i].dualHist[3*j] << ":";
			cout << " " << genes[i].dualHist[3*j+1] << ":";
			cout << " " << genes[i].dualHist[3*j+2] << " ";
			cout << " " << genes[i].subMTX[3*j] << ":";
			cout << " " << genes[i].subMTX[3*j+1] << ":";
			cout << " " << genes[i].subMTX[3*j+2] << endl;
		}
	}
	cout << "----------------------" << endl;
}


// information loaders
vector<Sequence> loadInfoFromInputFile(const char* infile, int N){
	FILE *inputFile = fopen(infile, "r");
	vector<Sequence> querySeqs;
	
	int maxLine = 5000;
	char line[maxLine];
	char delim = '\t';
	string oldQuery = "";
	string oldGene = "";
		
	while(fgets(line, maxLine, inputFile) != NULL){
		vector<string> elems = split(string(line), delim);
		string queryName = elems[12];
		string geneName = elems[13];
		IDnum geneGI = atoi(elems[14].c_str());
		float identity = atof(elems[2].c_str());
		float bitscore = atof(elems[11].c_str());
		
		if(oldQuery.compare(queryName) != 0){
			oldQuery.assign(queryName);
			Sequence querySequence;
			querySequence.seqName.assign(queryName);
			querySeqs.push_back(querySequence);
		}
		
		if(oldGene.compare(geneName) != 0){
			oldGene.assign(geneName);
			Gene queryGene;
			querySeqs.back().genes.push_back(queryGene);
		}
		
		if(querySeqs.back().genes.back().gis.size() < (unsigned int) N){
			querySeqs.back().genes.back().gis.push_back(geneGI);
			querySeqs.back().genes.back().identity.push_back(identity);
			querySeqs.back().genes.back().bitscore.push_back(bitscore);
		}
	}
	
	fclose(inputFile);

	return querySeqs;
}


// load gi->taxonID mapping information
void loadGI2TaxonLibFromFile(const char* gi2taxonFile, vector<Sequence> &QuerySeq){
	map<IDnum, IDnum> giHits;
	map<IDnum, IDnum>::iterator it;
	
	for(unsigned int index = 0; index < QuerySeq.size(); index++){
		for(unsigned int i = 0; i < QuerySeq[index].genes.size(); i++){
			for(unsigned int j = 0; j < QuerySeq[index].genes[i].gis.size(); j++){
				it = giHits.begin();
				giHits.insert(it, pair<IDnum, IDnum> (QuerySeq[index].genes[i].gis[j], 0));				
			}
		}
	}
	
	FILE *libFile = fopen(gi2taxonFile, "r");
	
	int maxLine = 1000;
	char line[maxLine];
	IDnum GI, taxonID;
	
	while(fgets(line, maxLine, libFile) != NULL){
		stringstream lineStream;
		lineStream << line;
		lineStream >> GI >> taxonID;
		if(giHits.count(GI) != 0){
			giHits.find(GI)->second = taxonID;
		}
	}
	fclose(libFile);
	
	for(unsigned int index = 0; index < QuerySeq.size(); index++){
		for(unsigned int i = 0; i < QuerySeq[index].genes.size(); i++){
			for(unsigned int j = 0; j < QuerySeq[index].genes[i].gis.size(); j++){
				IDnum GI = QuerySeq[index].genes[i].gis[j];
				if(GI > 0){
					QuerySeq[index].genes[i].taxonIDs.push_back(giHits.find(GI)->second);
				}else{
					QuerySeq[index].genes[i].taxonIDs.push_back(0);
				}
			}
		}
	}
}

float getHistPara(string dh, float identity){
	char delim = '\t';
	vector<string> elems = split(dh, delim);
	unsigned int index = 1000 - int(identity*10);
	try{
		return atof(elems[index].c_str());
	}catch(exception& e){
		return	1.0;
	}
}

// load gi->clstr mapping information
void loadGI2ClstrLibFromFile(const char* gi2clstrFile, vector<Sequence> &QuerySeq){
	map<IDnum, IDnum> gi2clstr;
	map<IDnum, IDnum>::iterator it;
	
	map<IDnum, vector<string> > dualHist;
	map<IDnum, vector<string> >::iterator dhit;
	
	map<IDnum, vector<float> > subMTX;
	map<IDnum, vector<float> >::iterator smit;
	
	gi2clstr.clear();
	dualHist.clear();
	subMTX.clear();

	for(unsigned int index = 0; index < QuerySeq.size(); index++){
		for(unsigned int i = 0; i < QuerySeq[index].genes.size(); i++){
			for(unsigned int j = 0; j < QuerySeq[index].genes[i].gis.size(); j++){
				it = gi2clstr.begin();
				
				gi2clstr.insert(it, pair<IDnum, IDnum> (QuerySeq[index].genes[i].gis[j], 0));	
			}			
		}
	}
	
	FILE *libFile = fopen(gi2clstrFile, "r");
	
	int maxLine = 10000;
	char line[maxLine];
	char delim = '\t';
	while(fgets(line, maxLine, libFile) != NULL){
		// first line, the clstr ID
		stringstream ss;
		ss << line;
		IDnum clstrID;
		int clstrSize;
		ss >> clstrID >> clstrSize;
		int numLines = clstrSize/10;
		if(clstrSize%10 != 0){
			numLines++;
		}
		
		//lines with all members of the gene cluster
		bool hasThisClstr = false;
		for(int lineNum = 0; lineNum < numLines; lineNum++){
			fgets(line, maxLine, libFile);
			vector<string> elems = split(string(line), delim);
			for(unsigned int index = 0; index < elems.size(); index++){
				IDnum GI = atoi(elems[index].c_str());
				if(gi2clstr.count(GI) != 0){
					hasThisClstr = true;
					gi2clstr.find(GI)->second = clstrID;
				}
			}
		}
		
		// handle line 3-6, phylum/genus/species parameters of this cluster
		if(hasThisClstr){
			dhit = dualHist.begin();
			smit = subMTX.begin();
			vector<string> dh;
			dh.clear();
			vector<float> sm;
			sm.clear();
			
			for(int i = 0; i < 3; i++){
				fgets(line, maxLine, libFile);
				dh.push_back(string(line));
			}
			dualHist.insert(dhit, pair<IDnum, vector<string> > (clstrID, dh));
			
			fgets(line, maxLine, libFile);
			vector<string> elems = split(string(line), delim);
			for(unsigned int j = 0; j < elems.size(); j++){
				try{
					float r = atof(elems[j].c_str());
					sm.push_back(r);
				}catch(exception& e){
					float r = -1.0;
					sm.push_back(r);
				}
			}
			subMTX.insert(smit, pair<IDnum, vector<float> > (clstrID, sm));
			
		}else{
			for(int i = 0; i < 4; i++){
				fgets(line, maxLine, libFile);
			}
		}
		
	}
	fclose(libFile);
	
	// load information onto QuerySeq
	
	try{
		for(unsigned int index = 0; index < QuerySeq.size(); index++){
			for(unsigned int i = 0; i < QuerySeq[index].genes.size(); i++){
				for(unsigned int j = 0; j < QuerySeq[index].genes[i].gis.size(); j++){
					IDnum GI = QuerySeq[index].genes[i].gis[j];
					IDnum clstrID = gi2clstr.find(GI)->second;
				
					if(clstrID == 0){   // in case the GI is not in lib;
						QuerySeq[index].genes[i].clusters.push_back(0);
						QuerySeq[index].genes[i].dualHist.push_back(-1);
						QuerySeq[index].genes[i].dualHist.push_back(-1);
						QuerySeq[index].genes[i].dualHist.push_back(-1);
						QuerySeq[index].genes[i].subMTX.push_back(-1);
						QuerySeq[index].genes[i].subMTX.push_back(-1);
						QuerySeq[index].genes[i].subMTX.push_back(-1);
						continue;
					}
				
					// regular case;
					vector<string> dh = dualHist.find(clstrID)->second;
					vector<float> sub = subMTX.find(clstrID)->second;
					
					float identity = QuerySeq[index].genes[i].identity[j];
					float phylumHist = getHistPara(dh[0], identity);
					float genusHist = getHistPara(dh[1], identity);
					float speciesHist = getHistPara(dh[2], identity);
					// load the parameters
					QuerySeq[index].genes[i].clusters.push_back(clstrID);
					QuerySeq[index].genes[i].dualHist.push_back(phylumHist);
					QuerySeq[index].genes[i].dualHist.push_back(genusHist);
					QuerySeq[index].genes[i].dualHist.push_back(speciesHist);
					
					for(int k = 0; k < 3; k++){
						QuerySeq[index].genes[i].subMTX.push_back(sub[k]);
					}
					
				}			
			}
		}
	}catch(exception& e){
		cout << "get an error: " << e.what() << endl;
	}
	//end of function
}

// add pertaining ranks taxonID to taxonPath of query sequences;
void addToSeqTaxonPaths(vector<IDRank> tPath, map<IDnum, PathNode*> &seqTaxonForest){
	PathNode *rootNode;
	map<IDnum, PathNode*>::iterator it;
	
	// if no root exists yet, then add a root to the forest;
	if(seqTaxonForest.size() == 0){
		rootNode = newPathNode();
		rootNode->category = 0;
		it = seqTaxonForest.begin();
		seqTaxonForest.insert(it, pair<IDnum, PathNode*> (0, rootNode));
	}else{
		rootNode = seqTaxonForest.find(0)->second;
	}
	// add taxonPath to the forest;
	PathNode *phylumNode = newPathNode();
	PathNode *genusNode = newPathNode();
	PathNode *speciesNode = newPathNode();
	
	for(unsigned int index = 0; index < tPath.size(); index++){
		IDnum taxonID = tPath[index].taxonID;
		string rank = tPath[index].rank;
		
		// discard the nodes that are not phylum/genus/species ranks;
		if((rank.compare("species") != 0) and (rank.compare("genus") != 0) and (rank.compare("phylum") !=0)){
			continue;
		}
	
		// either load the node onto the forest or point to pre-existing node;
		if(rank.compare("species") == 0){
			if(seqTaxonForest.count(taxonID) == 0){
				it = seqTaxonForest.begin();
				speciesNode->taxonID = taxonID;
				speciesNode->category = 3;
				speciesNode->likelihood = 0.0;
				seqTaxonForest.insert(it, pair<IDnum, PathNode*> (taxonID, speciesNode));
			}else{
				speciesNode = seqTaxonForest.find(taxonID)->second;
			}
		}else if(rank.compare("genus") == 0){
			if(seqTaxonForest.count(taxonID) == 0){
				it = seqTaxonForest.begin();
				genusNode->taxonID = taxonID;
				genusNode->category = 2;
				genusNode->likelihood = 0.0;
				seqTaxonForest.insert(it, pair<IDnum, PathNode*> (taxonID, genusNode));
			}else{
				genusNode = seqTaxonForest.find(taxonID)->second;
			}
		}else if(rank.compare("phylum") == 0){
			if(seqTaxonForest.count(taxonID) == 0){
				it = seqTaxonForest.begin();
				phylumNode->taxonID = taxonID;
				phylumNode->category = 1;
				phylumNode->likelihood = 0.0;
				seqTaxonForest.insert(it, pair<IDnum, PathNode*> (taxonID, phylumNode));
			}else{
				phylumNode = seqTaxonForest.find(taxonID)->second;
			}
		}
	} // end of for loop;
	
	// connect the nodes, hook to rootNode;
	if(phylumNode->prevNode == NULL){
		phylumNode->prevNode = rootNode;
	}
	if(genusNode->prevNode == NULL){
		genusNode->prevNode = phylumNode;
	}
	if(speciesNode->prevNode == NULL){
		speciesNode->prevNode = genusNode;
	}
	//end of function;
}

vector<IDnum> getTaxonIDAtThreeRanks(TaxonTree *tTree, IDnum leafTaxonID){
	vector<IDnum> taxonIDs;
	vector<IDRank> idr = taxonomyPathIDRank(tTree, leafTaxonID);
	IDnum phylum=0, genus=0, species=0;
	
	for(unsigned int index = 0; index < idr.size(); index++){
		IDnum taxonID = idr[index].taxonID;
		string rank = idr[index].rank;
		if(rank.compare("phylum") == 0){
			phylum = taxonID;
		}else if(rank.compare("genus") == 0){
			genus = taxonID;
		}else if(rank.compare("species") == 0){
			species = taxonID;
		}else{
			continue;
		}
	}
	
	taxonIDs.push_back(phylum);
	taxonIDs.push_back(genus);
	taxonIDs.push_back(species);
	
	return taxonIDs;
}

// calculate the likelihood of taxonomy for query sequences;
void likelihoodCal(TaxonTree *tTree, vector<Sequence> &QuerySeq){
	for(unsigned int seqIndex = 0; seqIndex < QuerySeq.size(); seqIndex++){
		// load all possible taxonomy paths onto query sequences;
		QuerySeq[seqIndex].seqTaxonForest.clear();
		
		for(vector<gene_st>::iterator git = QuerySeq[seqIndex].genes.begin();
				git != QuerySeq[seqIndex].genes.end(); ++ git){
			for(vector<IDnum>::iterator tit = git->taxonIDs.begin();
				tit != git->taxonIDs.end(); ++ tit){
				vector<IDRank> tPath = taxonomyPathIDRank(tTree, *tit);
				addToSeqTaxonPaths(tPath, QuerySeq[seqIndex].seqTaxonForest);
			}
				
		}
		
		// iterate through all matches and calculate the likelihoods;
		// extract pointers to nodes at different ranks, put them in vectors;
		vector<PathNode*> phylumNodes;
		vector<PathNode*> genusNodes;
		vector<PathNode*> speciesNodes;
		map<IDnum, PathNode*>::iterator forestIt;
		
		for(forestIt = QuerySeq[seqIndex].seqTaxonForest.begin(); 
			forestIt != QuerySeq[seqIndex].seqTaxonForest.end(); forestIt++){
			PathNode* node = forestIt->second;
			if(node->category == 0){
				continue;
			}else if(node->category == 1){
				phylumNodes.push_back(node);
			}else if(node->category == 2){
				genusNodes.push_back(node);
			}else if(node->category == 3){
				speciesNodes.push_back(node);
			}
		}
		
		//iterate through matches, and add up the scores;
		for(unsigned int geneIndex = 0; geneIndex < QuerySeq[seqIndex].genes.size(); geneIndex++){
			for(unsigned int i = 0; i < QuerySeq[seqIndex].genes[geneIndex].gis.size(); i++){
				
				IDnum leafTaxonID = QuerySeq[seqIndex].genes[geneIndex].taxonIDs[i];
				if(leafTaxonID == 0){
					continue;
				}
				
				vector<IDnum> taxonIDs = getTaxonIDAtThreeRanks(tTree, leafTaxonID);
				
				if(taxonIDs[0] == 0 or taxonIDs[1] ==0 or taxonIDs[2] == 0){
					continue;
				}
				
				PathNode* phylumNode = QuerySeq[seqIndex].seqTaxonForest.find(taxonIDs[0])->second;
				PathNode* genusNode = QuerySeq[seqIndex].seqTaxonForest.find(taxonIDs[1])->second;
				PathNode* speciesNode = QuerySeq[seqIndex].seqTaxonForest.find(taxonIDs[2])->second;
				
				float dhPhylum, dhGenus, dhSpecies;
				float smPhylum, smGenus, smSpecies;
				
				dhPhylum = QuerySeq[seqIndex].genes[geneIndex].dualHist[3*i];
				smPhylum = QuerySeq[seqIndex].genes[geneIndex].subMTX[3*i];
				
				dhGenus = QuerySeq[seqIndex].genes[geneIndex].dualHist[3*i+1];
				smGenus = QuerySeq[seqIndex].genes[geneIndex].subMTX[3*i+1];
				
				dhSpecies = QuerySeq[seqIndex].genes[geneIndex].dualHist[3*i+2];
				smSpecies = QuerySeq[seqIndex].genes[geneIndex].subMTX[3*i+2];
				
				phylumNode->likelihood += W10*dhPhylum + W20*smPhylum;
				genusNode->likelihood += W11*dhGenus + W21*smGenus;
				speciesNode->likelihood += W12*dhSpecies + W22*smSpecies;
			}
		}
		
		//iterate through nodes at different ranks in the taxonomy forest, and normalize scores into likelihoods;
		// sum of scores
		float phylumSum, genusSum, speciesSum;
		phylumSum = 0;
		genusSum = 0;
		speciesSum = 0;
		// phylum level;
		for(unsigned int index = 0; index < phylumNodes.size(); index++){
			phylumSum += phylumNodes[index]->likelihood;
		}
		
		for(unsigned int index = 0; index < phylumNodes.size(); index++){
			if(phylumSum != 0){
				phylumNodes[index]->likelihood /= phylumSum;
			}else{
				phylumNodes[index]->likelihood = 1.0;
			}
		}
		
		// genus level;
		for(unsigned int index = 0; index < genusNodes.size(); index++){
			genusSum += genusNodes[index]->likelihood;
		}
		
		for(unsigned int index = 0; index < genusNodes.size(); index++){
			if(genusSum != 0){
				genusNodes[index]->likelihood /= genusSum;
			}else{
				genusNodes[index]->likelihood = 1.0;
			}
		}
		
		//species level;
		for(unsigned int index = 0; index < speciesNodes.size(); index++){
			speciesSum += speciesNodes[index]->likelihood;
		}
		
		for(unsigned int index = 0; index < speciesNodes.size(); index++){
			if(speciesSum != 0){
				speciesNodes[index]->likelihood /= speciesSum;
			}else{
				speciesNodes[index]->likelihood = 1.0;
			}
		}
		
	}//end of for loop;
	
	//end of function;
}

// output results
void writeResultsToOutputFile(const char* outfile, TaxonTree *tTree, TaxonName *tName, 
									vector<Sequence> &QuerySeq, float thr){
	
	ofstream outputFile;
	outputFile.open(outfile, ios::out);
	
	for(unsigned int seqIndex = 0; seqIndex < QuerySeq.size(); seqIndex++){
		// start with species level;
		vector<PathNode*> speciesNodes;
		map<IDnum, PathNode*>::iterator forestIt;
		
		for(forestIt = QuerySeq[seqIndex].seqTaxonForest.begin(); 
			forestIt != QuerySeq[seqIndex].seqTaxonForest.end(); forestIt++){
			PathNode* node = forestIt->second;
			if(node->category == 3){
				speciesNodes.push_back(node);
			}
		}
		
		float maxSpeciesLLH = 0;
		IDnum speciesID;
		for(unsigned int index = 0; index < speciesNodes.size(); index++){
			float likelihood = speciesNodes[index]->likelihood;
			IDnum taxonID = speciesNodes[index]->taxonID;
			if(likelihood > maxSpeciesLLH){
				maxSpeciesLLH = likelihood;
				speciesID = taxonID;
			}
		}
		
		if(maxSpeciesLLH > thr){
			vector<NameRank> path = taxonomyPath(tTree, tName, speciesID);
			string pathString = taxonomyPathString(path);
			// write to file;
			outputFile << QuerySeq[seqIndex].seqName << "\tSpecies\t" << maxSpeciesLLH << "\t" << speciesID << endl;
			outputFile << pathString << endl;
			continue;
		}
		
		// no species level satisfies threshold, continue to genus level;
		vector<PathNode*> genusNodes;
		
		for(forestIt = QuerySeq[seqIndex].seqTaxonForest.begin(); 
			forestIt != QuerySeq[seqIndex].seqTaxonForest.end(); forestIt++){
			PathNode* node = forestIt->second;
			if(node->category == 2){
				genusNodes.push_back(node);
			}
		}
		
		float maxGenusLLH = 0;
		IDnum genusID;
		for(unsigned int index = 0; index < genusNodes.size(); index++){
			float likelihood = genusNodes[index]->likelihood;
			IDnum taxonID = genusNodes[index]->taxonID;
			if(likelihood > maxGenusLLH){
				maxGenusLLH = likelihood;
				genusID = taxonID;
			}
		}
		
		if(maxGenusLLH > thr){
			vector<NameRank> path = taxonomyPath(tTree, tName, genusID);
			string pathString = taxonomyPathString(path);
			// write to file;
			outputFile << QuerySeq[seqIndex].seqName << "\tGenus\t" << maxGenusLLH << "\t" << genusID << endl;
			outputFile << pathString << endl;
			continue;
		}
		
		// no genus level satisfies threshold, continue to phylum level;
		vector<PathNode*> phylumNodes;
		
		for(forestIt = QuerySeq[seqIndex].seqTaxonForest.begin(); 
			forestIt != QuerySeq[seqIndex].seqTaxonForest.end(); forestIt++){
			PathNode* node = forestIt->second;
			if(node->category == 1){
				phylumNodes.push_back(node);
			}
		}
		
		float maxPhylumLLH = 0;
		IDnum phylumID;
		for(unsigned int index = 0; index < phylumNodes.size(); index++){
			float likelihood = phylumNodes[index]->likelihood;
			IDnum taxonID = phylumNodes[index]->taxonID;
			if(likelihood > maxPhylumLLH){
				maxPhylumLLH = likelihood;
				phylumID = taxonID;
			}
		}
		
		if(maxPhylumLLH > thr){
			vector<NameRank> path = taxonomyPath(tTree, tName, phylumID);
			string pathString = taxonomyPathString(path);
			// write to file;
			
			outputFile << QuerySeq[seqIndex].seqName << "\tPhylum\t" << maxPhylumLLH << "\t" << phylumID << endl;
			outputFile << pathString << endl;
			continue;
		}
		
		// no phylum level satisfies threshold, mark as novel;
		// write to file as novel;
		
		outputFile << QuerySeq[seqIndex].seqName << "\tUnknown\tNA\tNA" << endl;
		outputFile << "NA" << endl;
	}
	
	outputFile.close();
}

