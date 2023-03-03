#pragma GCC optimize(3,"Ofast","inline")

#ifndef _SPECIES_H
#define _SPECIES_H

#include "utility.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class Species {

public:
	Species():nrows(0), ncols(0) {};
	Species(const char * infilename);
	MatrixXf Get_Abds();
	vector<string> Get_Sample_Names();
	vector<string> Get_Species_Names();
	friend ostream & operator << (ostream& os, const Species& M);
private:
	MatrixXf abds;
	int nrows, ncols;
	vector<string> samples;	
	vector<string> species;
};

Species::Species(const char * infilename) {
	ifstream infile(infilename, ifstream::in);
	if(!infile) {
		cerr << "Error: Cannot open file : " << infilename << endl;
		exit(0);
	}
	
	string buffer;
	getline(infile, buffer);
	stringstream strin(buffer);
	
	string sample_name;
	strin >> sample_name;
	while(strin >> sample_name) 
		samples.push_back(sample_name);
	ncols = (int)samples.size();

	vector<float> temp;
        string sp;
        float value;
        while(getline(infile, buffer)){
                stringstream strin(buffer);
                strin >> sp;
                species.push_back(sp);
                for(int j = 0; j < ncols; ++j) {
                        strin >> value;
                        temp.push_back(value);
                }
        }

	nrows = (int) species.size();
        abds.resize(nrows, ncols);

        for(int i = 0; i < nrows; ++i) {
                for(int j = 0; j < ncols; ++j) {
                	abds(i, j) = temp[i*ncols+j];
                }
        }
}

MatrixXf Species::Get_Abds() {
	return abds;
}

vector<string> Species::Get_Sample_Names() {
	return samples;
}

vector<string> Species::Get_Species_Names() {
        return species;
}

ostream & operator << (ostream& os, const Species& M) {
	for(int i = 0; i < M.nrows; ++i) {
		for(int j = 0; j < M.ncols; ++j) {
			os << M.abds(i, j) << "\t";
		}
		os << endl;
	}
	return os;
}

#endif
