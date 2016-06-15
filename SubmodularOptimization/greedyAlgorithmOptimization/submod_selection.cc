/*
 * Feature based submodular selection
 *
 *  Created on: May 26, 2016
 *      Author: Kai Wei (kaiwei@uw.edu)
 *
 *	The code implements the feature based submodular selection for selecting K-mers
 *
 *	Speed up with accelerated greedy algorithm (Minoux 1976)
 *  Objective g(S) = \sum_{f \in F} w_f g(m_f(S)), where F is the set of K-mers, m_f(S) defines the count of the K-mers f occurring in a set regions S, and w_f defines the weight for this K-mers.
 *
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"

#include "arguments.h"
#include "error.h"

#include <queue>

using namespace std;

#define SMALLEST_NUMBER -1e70
#define LARGEST_NUMBER 1e70

char* strInputMat=NULL;
char* strWeights=NULL;
char* strReference=NULL;
double per = 0;
int nUtt = 0;
char* strOutput=NULL;
char* costFile=NULL;
char* help=NULL;
int verb=0;
double alpha=1; // always set to 1
bool isList=true; // always set to true
int nFeatures = 400000; // might need to be modified if we deal with different data
double approx_greedy=0.9; // controls the relaxation of the greedy procedure. The smaller, the less of approximation to greedy, but the faster the code runs.  

struct Utterance{
     long int index;
     int num_uniq_wrds;
     int* digitWrds;
     float* featureVec;
     int tot_num_wrds;
};


Arg Arg::Args[]={

    Arg("graph", Arg::Req, strInputMat, "the graph file in ascii format",Arg::SINGLE),
    Arg("weight", Arg::Opt, strWeights, "the file defining the weight for each feature in ascii format",Arg::SINGLE),
    Arg("reference", Arg::Opt, strReference, "the file defining the reference for each genome region in ascii format",Arg::SINGLE),
    Arg("n", Arg::Req, nUtt, "the total number of regions",Arg::SINGLE),
    Arg("per", Arg::Req, per, "the percentage of data to be selected",Arg::SINGLE),
    Arg("out", Arg::Req, strOutput, "the output file",Arg::SINGLE),
    Arg("cost", Arg::Opt, costFile, "the cost file", Arg::SINGLE),
    Arg("alpha",Arg::Opt, alpha, "regularizer of cost", Arg::SINGLE),
    Arg("verb", Arg::Opt, verb, "verbosity",Arg::SINGLE),
    Arg("help", Arg::Help, help, "Print this message"),
    Arg()
};


class Increment{
	public: 
		Increment(){};
		Increment(double x, int i){ value=x; index=i; }
		bool operator< (const Increment&) const;
	
		int get_index() const { return index; }
		double get_value() const{ return value; }
	private:
		int index;
		double value;
};

bool Increment::operator< (const Increment& right) const
{
	return value < right.value;
}


long GetFileSize(const char * filename )
{
    struct stat statebuf;

    if ( stat( filename, &statebuf ) == -1 )
        return -1L;
    else
        return statebuf.st_size;
}


double concave_function(double K){
    //return sqrt(K); 
    return log(1+K); // use the log function as the concave function
}

double weighted_log_modular(struct Utterance* file, int nV, double* preCompute, int newi, int nFeatures, double intermediate_sum, vector<double> weightsFeature){
    double sum = 0;
    double diff;
    double temp;
    int num_wrds = file[newi].num_uniq_wrds;
    sum = intermediate_sum;
    double check = 0;
    for (int i =0; i<num_wrds; i++){
        temp = preCompute[file[newi].digitWrds[i]] + file[newi].featureVec[i];
        diff = weightsFeature[file[newi].digitWrds[i]]*(concave_function(temp) - concave_function(preCompute[file[newi].digitWrds[i]]));
        sum += diff;
    }
    return (sum);
    
}


double log_modular(struct Utterance* file, int nV, double* preCompute, int newi, int nFeatures, double intermediate_sum){
// this function implements the function evaluation of the form \sum_{f \in F} g (m_f(S)), where g is concave. 

    double sum = 0;
    double diff;
    double temp;
    int num_wrds = file[newi].num_uniq_wrds;
    sum = intermediate_sum;
    double check = 0;
    for (int i =0; i<num_wrds; i++){
        temp = preCompute[file[newi].digitWrds[i]] + file[newi].featureVec[i];
        diff = concave_function(temp) - concave_function(preCompute[file[newi].digitWrds[i]]);
        sum += diff;
    }
    return (sum);
}




int 
line2words(char *s, struct Utterance * utterance, int* feature_tot_cnt){ // process the input feature file
    int digitwrd[16000];
    float featureval[16000];
    int pos = 0;
    int a = 0;
    while (sscanf(s,"%d %f %n",&digitwrd[a], &featureval[a], &pos) == 2) {
         s += pos;
         a++;
    }
    for (int i=0; i<a; i++){
        if (digitwrd[i] >= nFeatures ){
            error("You need to enlarge the nFeatures, the index of the feature exceeds the total number of features to process\n");
        }
        feature_tot_cnt[digitwrd[i]] += featureval[i];
    }
    utterance->digitWrds = new int[a];
    utterance->featureVec = new float[a];
    memcpy(utterance->digitWrds, digitwrd, sizeof(int)*a);
    memcpy(utterance->featureVec, featureval, sizeof(float)*a);
    return a;
}

int main(int argc, char** argv){
	
	bool parse_was_ok = Arg::parse(argc,(char**)argv);
        if(!parse_was_ok){
          Arg::usage(); exit(-1);
        }
	
	int i,j;
	double tmpd; float tmpf;
	struct Utterance * file = new struct Utterance[nUtt]; //file is array of structure named utterance
    
	// read the list of distributed graphs
	vector<string> graphList;
	string tmp;
	ifstream iFile;
    	int* feature_tot_cnt = new int[nFeatures](); 
	if (isList){
		FILE *fp = NULL;
		char line2[300000]; //stores information in each line of input
		printf("Reading list of graph files from %s...\n", strInputMat);
		if ((fp = fopen(strInputMat, "rt")) == NULL)
		{
			printf("Error: Cannot open file %s", strInputMat);
			exit(-1);
		}
		long int lineno = 0;
		while ( fgets(line2,sizeof(line2),fp) != NULL){
			file[lineno].index = lineno;
			file[lineno].num_uniq_wrds = line2words(line2, &file[lineno], feature_tot_cnt);//line2words transforms input with fmt digwords:featurevals into initialization of structure "Utterance"         
			//printf("There are %d uniq triphones in %d utt\n",file[lineno].num_uniq_wrds, lineno);
			lineno ++;
			//printf("read in line %d\n", lineno);
		}
		fclose(fp);
		if (nUtt != lineno){
			printf("number of lines in input graph: %d\n", lineno);
			printf("error: number of points doesn't match with input number of points\n");
			exit(-1);
		}
	}
 
    // read the reference file 
    vector<string> refVec;
    if (strReference==NULL){
        for(int idx=0;idx<nUtt;idx++){
            refVec.push_back("");
        }
    }
    else{
        ifstream file(strReference);
        std::string s;
        while (std::getline(file, s))
        {
            refVec.push_back(s);
        }
    }

    // read the weight FILE
    vector<double> weightVec(nFeatures,1);
    if (strWeights==NULL){
        for (int idx=0;idx<nFeatures;idx++){
            weightVec[idx] = 1;
        }
    }
    else{
        // if the weight file has been defined, then try to read the file.
        // TODO: need to implement the method for reading the weight file.  
    }
    
	
	// read the list of costs
	vector<double> costList;
    	int *cost_number_wrds = new int[nUtt];
	double totalCost = 0;
	
	printf("Reading list of costs from %s...\n", costFile);
	iFile.open(costFile, ios::in);
	if (!iFile.is_open()){ // if the list of costs is not defined, use the uniform cost
		printf("Use unit cost for each utterance\n");
		for (i=0; i<nUtt; i++) {	
			costList.push_back(1.0);
			totalCost += 1.0;
		}	
	}
	else {
		for (i=0; i<nUtt; i++) {
			iFile >> tmpd;
			totalCost += tmpd;
			costList.push_back(tmpd);
		}
	}
	iFile.close();
	
	double maxCost = totalCost*(per/100);
	printf("Total cost of the dataset is %f.\n", totalCost);
	printf("Maximum budget for the subset is %f.\n", maxCost);
	
	int nV = nUtt ;
	int nRow;
	int count = 0; // row count
	FILE* fp;
	
    	int nSelected = 0 ; // number of utterances selected
	int* SelectedSet = new int[nV]; // array to store the selected utt ids
	bool* hashSelected = new bool[nV]; // an index table to show whether an utt is selected
	for(i=0;i<nV;i++) hashSelected[i]=false;	
    
	double* preCompute = new double[nFeatures]; //for speed-up, compute rho for each feature  = sqrt(sum over this feature values in all utts in currently selected subset)	
    	double intermediate_sum = 0;
    	for (i = 0; i < nFeatures; i++) { preCompute[i] = 0; }
    	double maxV;
    	double tmpV;
    	double preV = 0;
    	double objV = 0;
	double currentCost = 0;
    	int maxI;
	
 	clock_t start = clock();	
	// accelerated greedy algorithm implementation
	printf("Start accelerated greedy algorithm\n");
	priority_queue <Increment> rho;
	// initilize the priority queue
	printf("Initilize priority queue\n");
	for (i = 0; i < nV; i++) {
		// evaluate every i \in ground_set V
		//rho.push(Increment(log_modular(file,nV,preCompute,i,nFeatures, intermediate_sum)/pow(costList[i],alpha),i)); // I think it makes more sense if we set alpha to be 1, earlier, we set it to be 0, which does not yield good result. It's necessary to also try out using number of words as normalizer in this function, which I think should make more sense. 
		rho.push(Increment(weighted_log_modular(file,nV,preCompute,i,nFeatures, intermediate_sum, weightVec)/pow(costList[i],alpha),i)); // I think it makes more sense if we set alpha to be 1, earlier, we set it to be 0, which does not yield good result. It's necessary to also try out using number of words as normalizer in this function, which I think should make more sense. 
        // BTW, currently, I think we should turn off modular trick
	}
	cout << endl;

	int sort_cnt = 0;
	
	while (! rho.empty()) {
		int topid = rho.top().get_index();
		rho.pop();
		//maxV = sum_sqrt_modular(file,nV,preCompute,topid,nFeatures,intermediate_sum);
		//maxV = log_modular(file,nV,preCompute,topid,nFeatures,intermediate_sum);
		maxV = weighted_log_modular(file,nV,preCompute,topid,nFeatures,intermediate_sum, weightVec);
		double newV = (maxV - preV)/pow(costList[topid],alpha);
		if (verb >= 5) printf("max gain = %.6e, rho->top() = %.6e\n",newV,rho.top().get_value());
		if (newV < rho.top().get_value()*approx_greedy) {
			rho.push(Increment(newV, topid)); // if condition not satisfied, push back and re-sort
			sort_cnt++;
			if (verb >= 10) printf("Condition not met, re-sorting queue (%d iterations)\n",sort_cnt);
		}
		else {
			// guaranteed to be optimal because of submodularity
		    	hashSelected[topid] = true;
            		SelectedSet[nSelected++] = topid;
			objV += newV;
			currentCost += costList[topid];
			printf("Selecting %dth sample (%d selected, current normalized increment = %.6e), curCost/budget = %.6e/%.6e, preV = %.10e, curV = %.10e, num_resort=%d\n", topid, nSelected, newV, currentCost, maxCost, preV, maxV, sort_cnt);
            
			
            		preV = maxV;
            		intermediate_sum = preV;
            
            		double diff_sum = 0;
            		for (i = 0; i < file[topid].num_uniq_wrds; i++){ // only compute updates for precompute that changes
                		preCompute[file[topid].digitWrds[i]] += file[topid].featureVec[i];
            		}
	            //intermediate_sum += diff_sum; // update function evaluation on currently selected subset
          		  /*
	        for (i = 0; i < nFeatures; i++) {
                preCompute[i] += fKernel[topid][i];
			}*/
			sort_cnt = 0;
			if (currentCost >= maxCost) 
				break;
		}    
	}
        clock_t end = clock();

	cout << "In total running: " << (double) (end - start) / CLOCKS_PER_SEC << " seconds" << endl;	
	//printf("Total number of resorting is %d\n", sort_cnt);
	FILE* ofp = fopen(strOutput,"w");
	for(i=0;i<nSelected;i++)
		fprintf(ofp,"%d ",SelectedSet[i]);	
	fclose(ofp);

    ofstream myfile;
    //myfile.open(strOutput+".ref");
    string outref("");
    i=0;
    while(strOutput[i]!=NULL){
        outref.push_back(strOutput[i]);
        i++;
    }
    outref = outref+".ref";
    myfile.open(outref.c_str());
    for (int idx=0;idx<nSelected;idx++){
        myfile << refVec[SelectedSet[idx]] << endl;
    }
    myfile.close();
	return 0;
	
}
