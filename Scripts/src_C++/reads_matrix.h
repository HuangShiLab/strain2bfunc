#include "utility.h"

class _Reads_Matrix {

public:
	
	_Reads_Matrix(){}

	_Reads_Matrix(const char * reads_file, int file_type){
		switch(file_type) {
			case 0:
				break;
			case 1:
				break;
			case 2:
				break;
			case 3:
				break;
		}
	}

	void Load_Fasta(const char * reads_file, int sample_number = 0);
	void Load_Fastq(const char * reads_file, int sample_number = 0);
	void Load_Fasta_List(const char * list_file);
	void Load_Fastq_List(const char * list_file);

private:
	map<string, int> reads_index;
	vector<string> reads_vector;
	vector<vector<int> > reads_matrix;
};

_Reads_Matrix::_Reads_Matrix(const char * infilename, int file_type){
	switch(file_type) {
		case 0:
			Load_Fasta(infilename);
			break;
		case 1:
			Load_Fastq(infilename);
			break;
		case 2:
			Load_Fasta_List(infilename);
			break;
		case 3:
			Load_Fastq_List(infilename);
			break;
	}
}

void _Reads_Matrix::Load_Fasta(const char * infilename, int sample_number){

	ifstream infile(infilename, ifstream::in);

	if (!infile){
		cerr << "Error: Cannot open file : " << infilename << endl;
		return 0;
	}

	string buffer;

	unsigned int count = 0;
	int index = 0;

	while (getline(infile, buffer)){

		if (buffer[0] != '+' && buffer.size() == 1) {
			
			map<string, int>::iterator iter;
			iter = reads_index.find(buffer);
			if(iter != reads_index.end()) {
				index = reads_index[buffer];
			}
			else {
				reads_index.insert(std::pair<string, int>(buffer, reads_vector.size()));
				reads_vector.push_back(buffer);
			}
		}
	}

	infile.close();
	infile.clear();

	return count;
}	
