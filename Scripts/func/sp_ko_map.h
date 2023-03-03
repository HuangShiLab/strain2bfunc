#pragma GCC optimize(3,"Ofast","inline")

#ifndef _SPKOMAP_H
#define _SPKOMAP_H

#include "utility.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class Sp_ko_map {

public:
	Sp_ko_map():nrows(0), ncols(0) {};
	Sp_ko_map(const char * infilename, const vector<string> & sa_species);
	MatrixXf Get_Map();
	vector<string> Get_KO_Names();
	vector<string> Get_Species_Names();
	friend MatrixXf operator * (Sp_ko_map& m, Species& s);
	friend ostream & operator << (ostream& os, const Sp_ko_map& m);
private:
	MatrixXf maps;
	int nrows, ncols;
	vector<string> kos;	
	vector<string> species;
};

Sp_ko_map::Sp_ko_map(const char * infilename, const vector<string> & sa_species) {
	ifstream infile(infilename, ifstream::in);
	if(!infile) {
		cerr << "Error: Cannot open file : " << infilename << endl;
		exit(0);
	}

	species = sa_species;
	map<string, int> sp_index;
        ncols = (int)species.size();
        for(int i = 0; i < ncols; ++i) {
                sp_index[species[i]] = i;
        }

	string buffer;
	getline(infile, buffer);
	stringstream strin(buffer);
	
	string ko_number;
        strin >> ko_number;
        while(strin >> ko_number)
                kos.push_back(ko_number);
        nrows = (int)kos.size();

	maps = MatrixXf::Zero(nrows, ncols);

	/*
        maps.resize(nrows, ncols);
	for(int i = 0; i < nrows; ++i) {
		for(int j = 0; j < ncols; ++j) {
			maps(i, j) = 0;
		}
	}
	*/	

	string sp_name;
        float abd;
        while(getline(infile, buffer)){ //读进来之后要进行转置，即行名为ko number，列名为species name
                stringstream strin(buffer);
                strin >> sp_name;
		map<string, int>::iterator index = sp_index.find(sp_name);
		if (index == sp_index.end())
			continue;
                for(int j = 0; j < nrows; ++j) {
                        strin >> abd;
                        maps(j, index->second) = abd;
                }
        }
}

MatrixXf Sp_ko_map::Get_Map() {
	return maps;
}

vector<string> Sp_ko_map::Get_KO_Names() {
	return kos;
}

vector<string> Sp_ko_map::Get_Species_Names() {
        return species;
}

MatrixXf operator * (Sp_ko_map& m, Species& s){
	MatrixXf t = m.maps * s.Get_Abds();
	return t;
}

ostream & operator << (ostream& os, const Sp_ko_map& m) {
	for(int i = 0; i < m.nrows; ++i) {
		for(int j = 0; j < m.ncols; ++j) {
			os << m.maps(i, j) << "\t";
		}
		os << endl;
	}
	return os;
}

#endif
