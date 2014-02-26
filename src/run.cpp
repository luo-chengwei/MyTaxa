/*

	This file is part of MyTaxa by Chengwei Luo (luo.chengwei@gatech.edu)
    Konstantinidis Lab, Georgia Institute of Technology, 2013

*/

#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <sys/stat.h>

#include "run.h"

using namespace std;

////////////////////////// CLASS & STRUCTS ////////////////////////
// prints usage of Arguseq
void printUsage()
{
	cout << "#############################################################################################" << endl;
	cout << "MyTaxa: an advanced taxonomy classifier for metagenomic and genomic sequences" <<endl;
	cout << "Version: " << VERSION_NUMBER << ".";
	cout << RELEASE_NUMBER << "." << UPDATE_NUMBER << endl;
	cout << "Usage:" << endl;
	cout << "MyTaxa <input file> <output file> ";
	cout << "<score cutoff> <num of hit to use>" << endl;
	cout << "## [Format of input file]:" << endl;
	cout << "\tBased on blast -m 8 output format, for each blast-like output line," << endl;
	cout << "\tadd additional 3 tab delimited columns to each line:" << endl;
	cout << "\t[Query sequence name] [Gene name] [protein GI number]" << endl;
	cout << "#############################################################################################" << endl;
}

//defines an exception class for parseCommandArgs;
class myException: public exception{
	virtual const char* what() const throw(){
    	return "Command line arguments error";
    }
} myex;

//class that defines the arguments Arguseq uses;
class commandArgs{
public:
	const char* inputFile;
	const char* outputFile;
	float scoreThr;
	int numHits;
		
	void printArgs(){
		cout << "## The input file is: " << inputFile << endl;
		cout << "## The output will be stored at: " << outputFile << endl;
		cout << "## The output score cutoff is: " << scoreThr <<endl;
		cout << "## The number of hits used per gene is: "<< numHits << endl;
	}
}Args;

void initArgs(int argc, char** argv, commandArgs &Args){
	if (argc != 5) {
		throw myex;
	}else{
		try{
			Args.inputFile = realpath(argv[1], NULL);
			if(Args.inputFile == NULL){
				throw myex;
			}
			Args.outputFile = argv[2];
			Args.scoreThr = atof(argv[3]);
			Args.numHits = atoi(argv[4]);
		}catch(exception &e){
			cerr << "Argument error: " << e.what() << endl;
			exit(1);
		}
	}
}

//class that defines all the files in the ./db directory
class databaseFiles{
public:
	const char* taxonTreeFile;
	const char* taxonSciNameFile;
	const char* geneTaxonFile;
	const char* geneInfoFile;
	
	void initDBFiles(char *progPath){
		string progString = string(progPath);
		unsigned int pos = progString.rfind("/");
		string dbPathString = progString.substr(0, pos) + "/db/";
		string taxonTreeFileString = dbPathString + "ncbiNodes.lib";
		string taxonSciNameFileString = dbPathString + "ncbiSciNames.lib";
		string geneTaxonFileString = dbPathString + "geneTaxon.lib";
		string geneInfoFileString = dbPathString + "geneInfo.lib";
		
		taxonTreeFile = realpath(taxonTreeFileString.c_str(), NULL);
		taxonSciNameFile = realpath(taxonSciNameFileString.c_str(), NULL);
		geneTaxonFile = realpath(geneTaxonFileString.c_str(), NULL);
		geneInfoFile = realpath(geneInfoFileString.c_str(), NULL);
	}
	
}dbFiles;


////////////////////////// MAIN ///////////////////////
int main(int argc, char** argv){
	//init the argument for the run
	try{
		initArgs(argc, argv, Args);
		Args.printArgs();
	}catch(exception& e){
		cout<< e.what()<<endl;  // if error happens, reports here;
		printUsage();
		return 1;	
	}
	
	//load all the ./db file vars;
	dbFiles.initDBFiles(argv[0]);
	
	//load ncbi taxonomy libs
	
	cout << "Loading NCBI taxonomy information..."<<endl;
	TaxonTree *tTree;
	TaxonName *sciName;
	tTree = importTaxonTreeFromFile(dbFiles.taxonTreeFile);
	sciName = importTaxonNameFromFile(dbFiles.taxonSciNameFile);
	cout << "Done!" << endl;
	
	//  read input file, load all gi# into vector<IDnum> gis, and initialize
	//  the vector<Sequence*> querySequences; 
	cout << "Loading input file..." << endl;
	vector<Sequence> QuerySeq;
	QuerySeq.clear();
	QuerySeq = loadInfoFromInputFile(Args.inputFile, Args.numHits);
	cout << "Done!" << endl;	

	// load pre-calculated parameters
	// step 1, load GI->taxonID
	cout << "Loading gi2taxonID library..." << endl;
	loadGI2TaxonLibFromFile(dbFiles.geneTaxonFile, QuerySeq);
	cout << "Done!" << endl;
	
	
	// step 2, load GI->gene cluster
	cout << "Loading gene cluster information and parameters..." << endl;
	loadGI2ClstrLibFromFile(dbFiles.geneInfoFile, QuerySeq);
	cout << "Done!" << endl;
	
	// step 3, calculate the taxonomy for each query sequence.
	cout << "Calculating likelihoods of taxonomy affiliations..." << endl;
	likelihoodCal(tTree, QuerySeq);
	cout << "Done!" << endl;
	
	// output results
	cout << "Outputting results..." << endl;
	writeResultsToOutputFile(Args.outputFile, tTree, sciName, QuerySeq, Args.scoreThr);
	cout << "Done!" << endl;
	
	// clean up;
	cout << "Cleaning up..." << endl;
	destroyTaxonTree(tTree);
	destroyTaxonName(sciName);
	
	cout << "All finished, results are stored in " << Args.outputFile << endl;
	cout << "Bye!" << endl;
	
	return 0;
}
