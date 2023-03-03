#include "species.h"
#include "sp_ko_map.h"

using namespace std;

void printhelp();
int Parse_Para(int argc, char * argv[]);
void write_file(MatrixXf m, vector<string>& ko_names, vector<string>& sample_names, const char* outfilename);

string sp_abd_file;
string sp_ko_map_file;
string ko_abd_file = "ko.abd";

int main(int argc, char *argv[]) {

	Parse_Para(argc, argv);

	Species species(sp_abd_file.c_str());
	vector<string> sp_names = species.Get_Species_Names();
	vector<string> sample_names = species.Get_Sample_Names();

        Sp_ko_map sp_ko_map(sp_ko_map_file.c_str(), sp_names);
	vector<string> ko_names = sp_ko_map.Get_KO_Names();

	MatrixXf ko_abds = sp_ko_map * species;

	write_file(ko_abds, ko_names, sample_names, ko_abd_file.c_str());

	return 0;
}

void printhelp(){
	cout << "Compute the KO abundance table" << endl;
	cout << "Usage: " << endl;
	cout << "calculate-ko-abd [Option] Value" << endl;
	cout << "Options:" << endl;
	cout << "[Input options, required]" << endl;
	cout << "\t-i Input species abundance file" << endl;
	cout << "\t-m Input map file of species and kos" << endl;
	cout << "[Output options]" << endl;
	cout << "\t-o Output kos abundance file, default is \"ko.abd\"" << endl;
	cout << "[Other options]" << endl;
	cout << "\t-h Help" << endl;
	exit(0);
}

int Parse_Para(int argc, char * argv[]){
	int i = 1;
	if (argc == 1) printhelp();
	while(i < argc){
		if (argv[i][0] != '-') {
			cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
			exit(0);
		};
		switch(argv[i][1]){
			case 'i': sp_abd_file = argv[i+1]; break;
			case 'm': sp_ko_map_file = argv[i+1]; break;
			case 'o': ko_abd_file = argv[i+1]; break;
			case 'h': printhelp(); break;
			default :
			cerr << "Error: Unrec argument " << argv[i] << endl;
				printhelp();
				break;
		}
		i += 2;
	}
	return 1;
}

void write_file(MatrixXf m, vector<string>& ko_names, vector<string>& sample_names, const char* outfilename) {
	ofstream outfile(outfilename, ofstream::out);
	if (!outfile){
		cerr << "Error: Cannot open output file: " << outfilename << endl;
		return;
	}
	
	int ncols = (int) sample_names.size();
	int nrows = (int) ko_names.size();
	outfile << "KO";
	for(int i = 0; i < ncols; ++i) {
		outfile << "\t" << sample_names[i];
	}
	outfile << endl;

	for(int i = 0; i < nrows; ++i) {
		outfile << ko_names[i];
		for(int j = 0; j < ncols; ++j) {
			outfile << "\t" << m(i, j);
		}
		outfile << endl;
	}	

	outfile.close();
	outfile.clear();
}
