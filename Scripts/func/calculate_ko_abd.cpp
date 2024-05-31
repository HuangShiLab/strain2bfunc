#include "species.h"
#include "sp_ko_map.h"

using namespace std;

void printhelp();
int Parse_Para(int argc, char * argv[]);
MatrixXf Count_to_Abd(const MatrixXf & count_mat, vector<string> & ko_names);
void write_file(const MatrixXf & m, vector<string>& ko_names, vector<string>& sample_names, const char* outfilename);

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

	MatrixXf ko_count = sp_ko_map * species;

	MatrixXf ko_abds = Count_to_Abd(ko_count, ko_names);

	write_file(ko_abds, ko_names, sample_names, ko_abd_file.c_str());

	return 0;
}

void printhelp(){
	cout << "Compute the KO abundance table" << endl;
	cout << "Usage: " << endl;
	cout << "calculate-ko-abd [Option] Value" << endl;
	cout << "Options:" << endl;
	cout << "[Input options, required]" << endl;
	cout << "\t-i Input strain abundance file" << endl;
	cout << "\t-m Input map file of strains and kos" << endl;
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

MatrixXf Count_to_Abd(const MatrixXf & count_mat, vector<string> & ko_names) {
    // Calculate the sum of each column
    VectorXf col_sum = count_mat.colwise().sum();

    // Avoid division by zero
    for (int i = 0; i < col_sum.size(); ++i) {
        if (col_sum(i) == 0) {
            // Set the column sum to a default value (e.g., 1)
            col_sum(i) = 1;
        }
    }

    // Convert to abundance
    MatrixXf abundance_mat = count_mat.array().rowwise() / col_sum.transpose().array();

    // Delete KO entries that are absent in all samples (i.e., rows with a sum of zero).
    // Calculate row sums
    VectorXf row_sums = abundance_mat.rowwise().sum();
    
    // Find rows with non-zero sums
    Array<bool, Dynamic, 1> nonzero_rows = (row_sums.array() != 0);
    int nonzero_rows_count = 0;
    for(int i = 0; i < nonzero_rows.size(); i++) {
	    if(nonzero_rows(i)) nonzero_rows_count++;
    }
    
    vector<string> new_ko_names;
    
    // Filter out rows with non-zero sums
    MatrixXf filtered_matrix(nonzero_rows_count, abundance_mat.cols());
    int filtered_row_idx = 0;
    for (int i = 0; i < nonzero_rows.size(); i++) {
        if (nonzero_rows(i)) {
            filtered_matrix.row(filtered_row_idx) = abundance_mat.row(i);
            new_ko_names.push_back(ko_names[i]);
	    filtered_row_idx += 1;
        }
    }

    ko_names = new_ko_names;

    return filtered_matrix;
}

void write_file(const MatrixXf & m, vector<string>& ko_names, vector<string>& sample_names, const char* outfilename) {
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
