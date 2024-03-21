// Updated at March 8, 2024
// Code by ZHANG Yufeng
// Faculty of dentistry, HKU

#include "pipeline.h"
#include "utility.h"
#include <omp.h>

using namespace std;

int main(int argc, char *argv[])
{

    Parse_Para(argc, argv);

    char command[BUFFER_SIZE];
    //Parallel for R
    vector<string> command_parallel_scripts;
    vector<string> command_parallel_titles;

    vector<string> command_parallel_bdiversity_scripts;
    vector<string> command_parallel_bdiversity_titles;

    Check_Path(Out_path.c_str(), 1);

    //open script file
    ofstream outscript((Out_path + "/scripts.sh").c_str(), ofstream::out);
    if (!outscript)
    {
        cerr << "Warning: Cannot open output file : " << Out_path << "/scripts.sh" << endl;
    }

    cout << "Strain2bFunc Version: " << Version << endl;
    outscript << "#Strain2bFunc Version: " << Version << endl;

    cout << "The reference sequence database is ";
    //cout << Database.Get_Description() << endl; // NEED to update: the versions of DBs used by strain2bfunc

    outscript << "#The reference sequence database is ";
    //outscript << Database.Get_Description() << endl; // NEED to update: the versions of DBs used by strain2bfunc

    //check metadata
    if (Load_ID(Meta_file.c_str(), Ids, 1) == 0)
    {
        string error_info = "Error: Please check the Meta data file (-m): at least contains 1 columns of Sample ID";
        cerr << error_info << endl;
        Echo_Error(error_info.c_str(), Error_file.c_str());
        return 0;
    }

    if (!Check_Ids(Ids))
    {
        string error_info = "Warning: Sample ID starts by number";
        cerr << error_info << endl;
        Echo_Error(error_info.c_str(), Error_file.c_str());
        //return 0;
    }

    Input_sam_num = Ids.size();
    vector<string> species_list;
    int species_count;
    string species;

    Check_Path(Temp_dir.c_str(), 1);

    cout << "mkdir -p " << Out_path << endl;
    sprintf(command, "mkdir -p %s", Out_path.c_str());
    outscript << command << endl;
    Run_With_Error(command, "Make directory", tmpError_file.c_str());

    switch (Step) {
        //Start from step 0
        case 0:
            //Step 0: 2bRAD-M
            
            if (Load_List(Seq_list_file.c_str(), Seq_files, List_prefix) == 0)
            {
                string error_info = "Error: Please check the sequence list file (-i) or the list path prefix (-p)";
                cerr << error_info << endl;
                Echo_Error(error_info.c_str(), Error_file.c_str());
                return 0;
            }
            
            //profiling
            cout << endl << "Microbial Community species-level profiling" << endl;
            outscript << endl << "#Microbial Community species-level profiling" << endl;
            
            cout << "perl " << path_2bRADM << "/bin/2bRADM_Pipline.pl -t " << format << " -l " << Seq_list_file << " -d " << database_path << "/2B-RAD-M-ref_db_GTDB -o " << Out_path <<  "/Species_results -qc no -gsc 5" << endl;
            sprintf(command, "perl %s/bin/2bRADM_Pipline.pl -t %d -l %s -d %s/2B-RAD-M-ref_db_GTDB -o %s/Species_results -qc no -gsc 5", path_2bRADM.c_str(), format, Seq_list_file.c_str(), database_path.c_str(), Out_path.c_str());
            outscript << command << endl;
            Run_With_Error(command, "2bRAD-M", tmpError_file.c_str());
            
            Taxa_list_file = Out_path + "/Species_results/list/BcgI.list";
            Table_file = Out_path + "/Species_results/quantitative/Abundance_Stat.all.xls";
            
            //Step 0 finished
            
        //Start from step 1
        case 1:
            //Step 1: Select species

            if (Taxa_list_file.size() == 0)
            {
                string error_info = "Error: Please check the taxa list (-l)";
                cerr << error_info << endl;
                Echo_Error(error_info.c_str(), Error_file.c_str());
                return 0;
            }
            
            cout << "Rscript " << path << "/Scripts/strain2b/make_species_list.R -i " << Table_file << " -t " << sp_thres << " -o " << Out_path << "/species_list.txt" << endl;
            sprintf(command, "Rscript %s/Scripts/strain2b/make_species_list.R -i %s -t %f -o %s/species_list.txt", path.c_str(), Table_file.c_str(), sp_thres, Out_path.c_str());
            outscript << command << endl;
            Run_With_Error(command, "make_species_list", Error_file.c_str());
            
            Species_list_file = Out_path + "/species_list.txt";
            //Step 1 finished

        //Start from step 2
        case 2:
            //Step 2: Strain-level profiling
            switch(Mode) {
                //Multiple-species merged
                case 1:
                    //Mkdir
                    cout << "mkdir -p " << Out_path << "/strain_results" << endl;
                    sprintf(command, "mkdir -p %s/strain_results", Out_path.c_str());
                    Run_With_Error(command, "Make directory of strain-level profiling results", Error_file.c_str());
                    outscript << command << endl; //We need to create the output directory before writing the scripts file

                    //Strain-level profiling
                    cout << "Rscript " << path << "/Scripts/strain2b/strain_pipeline.R -l " << Taxa_list_file << " -s " << Species_list_file << " -m 0 -d " << database_path << "/copy_number_matrix_0.001 -o " << Out_path << "/strain_results" << endl;
                    sprintf(command, "Rscript %s/Scripts/strain2b/strain_pipeline.R -l %s -s %s -m 0 -d %s/copy_number_matrix_0.001 -o %s/strain_results", path.c_str(), Taxa_list_file.c_str(), Species_list_file.c_str(), database_path.c_str(), Out_path.c_str());
                    outscript << command << endl;
                    Run_With_Error(command, "Strain-level profiling", Error_file.c_str());
                    
                    break;
                    
                //Multiple-species separated (Mode == 0)
                default:
                    species_count = Load_List(Species_list_file.c_str(), species_list);
                    
                    if(species_count == 0) {
                        string error_info = "There is no suitable species for strain-level profiling";
                        cerr << error_info << endl;
                        Echo_Error(error_info.c_str(), Error_file.c_str());
                        return 0;
                    }
                    
                    omp_set_num_threads(Min(species_count, Coren));
                    #pragma omp parallel for private(species, command)
                    for(int i = 0; i < species_count; i++) {
                        species = species_list[i];
                        
                        //Mkdir
                        cout << "mkdir -p " << Out_path << "/strain_results/" << species << endl;
                        sprintf(command, "mkdir -p %s/strain_results/%s", Out_path.c_str(), species.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Make directory of strain-level profiling results", Error_file.c_str());
                        
                        //Generate species list for each selected species
                        cout << "echo " << species << " > " << Out_path.c_str() << "/strain_results/" << species << "_list.txt" << endl;
                        sprintf(command, "echo %s > %s/strain_results/%s_list.txt", species.c_str(), Out_path.c_str(), species.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Generate Species list for strain-level profiling", Error_file.c_str());
                        
                        //Strain-level profiling
                        cout << "Rscript " << path << "/Scripts/strain2b/strain_pipeline.R -l " << Taxa_list_file << " -s " << Out_path << "/strain_results/" << species << "_list.txt -m 0 -d " << database_path << "/copy_number_matrix_0.001 -o " << Out_path << "/strain_results/" << species << endl;
                        sprintf(command, "Rscript %s/Scripts/strain2b/strain_pipeline.R -l %s -s %s/strain_results/%s_list.txt -m 0 -d %s/copy_number_matrix_0.001 -o %s/strain_results/%s", path.c_str(), Taxa_list_file.c_str(), Out_path.c_str(), species.c_str(), database_path.c_str(), Out_path.c_str(), species.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "strain-level profiling", Error_file.c_str());
                    }
                    break;
            }
            //Step 2 finished
            
            //Step 3: Function Prediction
            if (Is_func) {
                
                cout << endl << "Function Prediction" << endl;
                outscript << endl << "#Function Prediction" << endl;
                
                switch(Mode) {
                    //Multiple-species merged
                    case 1:
                        //Mkdir
                        cout << "mkdir -p " << Out_path << "/ko_results" << endl;
                        sprintf(command, "mkdir -p %s/ko_results", Out_path.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Make directory of function profiling results", Error_file.c_str());
                        
                        //Function profiling
                        cout << path << "/Scripts/func/calculate_ko_abd -i " << Out_path << "/strain_results/strain_level_abd.txt -m " << database_path << "/genome_to_ko.tsv -o " << Out_path << "/ko_results/ko_abd.txt" << endl;
                        sprintf(command, "%s/Scripts/func/calculate_ko_abd -i %s/strain_results/strain_level_abd.txt -m %s/genome_to_ko.tsv -o %s/ko_results/ko_abd.txt", path.c_str(), Out_path.c_str(), database_path.c_str(), Out_path.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Function profiling", Error_file.c_str());
                        
                        break;
                        
                    //Multiple-species separated (Mode == 0)
                    default:
                        if(species_count == 0) {
                            string error_info = "There is no suitable species for function profiling";
                            cerr << error_info << endl;
                            Echo_Error(error_info.c_str(), Error_file.c_str());
                            return 0;
                        }
                        
                        omp_set_num_threads(Min(species_count, Coren));
                        #pragma omp parallel for private(species, command)
                        for(int i = 0; i < species_count; i++) {
                            species = species_list[i];
                            
                            //Mkdir
                            cout << "mkdir -p " << Out_path << "/ko_results/" << species << endl;
                            sprintf(command, "mkdir -p %s/ko_results/%s", Out_path.c_str(), species.c_str());
                            outscript << command << endl;
                            Run_With_Error(command, "Make directory of function profiling results", Error_file.c_str());
                            
                            //Function profiling for each selected species
                            cout << path << "/Scripts/func/calculate_ko_abd -i " << Out_path << "/strain_results/" << species << "/strain_level_abd.txt -m " << database_path << "/genome_to_ko.tsv -o " << Out_path << "/ko_results/" << species << "/ko_abd.txt" << endl;
                            sprintf(command, "%s/Scripts/func/calculate_ko_abd -i %s/strain_results/%s/strain_level_abd.txt -m %s/genome_to_ko.tsv -o %s/ko_results/%s/ko_abd.txt", path.c_str(), Out_path.c_str(), species.c_str(), database_path.c_str(), Out_path.c_str(), species.c_str());
                            outscript << command << endl;
                            Run_With_Error(command, "Function profiling", Error_file.c_str());
                        }
                        break;
                }
                
            }
            //Step 3 finished
            
            //Step 4: Data analysis
            cout << endl << "Data Analysis" << endl;
            outscript << endl << "#Data Analysis" << endl;
            
            switch(Mode) {
                //Multiple-species merged
                case 1:
                    //prefix of analysis-results-file name
                    //prefix_name = "default";
                    
                    //Strain data analysis
                    //Mkdir
                    cout << "mkdir -p " << Out_path << "/strain_data_analysis_results" << endl;
                    sprintf(command, "mkdir -p %s/strain_data_analysis_results", Out_path.c_str());
                    outscript << command << endl;
                    Run_With_Error(command, "Make directory of strain data analysis results", Error_file.c_str());
                    
                    //Data analysis
                    cout << "sh " << path << "/Scripts/analysis/data_analysis.sh " << Out_path << "/strain_results/strain_level_abd.txt " << Meta_file << " " << Out_path << "/strain_data_analysis_results/ dist.txt " << prefix_name << endl;
                    sprintf(command, "sh %s/Scripts/analysis/data_analysis.sh %s/strain_results/strain_level_abd.txt %s %s/strain_data_analysis_results dist.txt %s", path.c_str(), Out_path.c_str(), Meta_file.c_str(), Out_path.c_str(), prefix_name.c_str());
                    outscript << command << endl;
                    Run_With_Error(command, "Strain data analysis", Error_file.c_str());
                    
                    //Function data analysis
                    if (Is_func) {
                        //Mkdir
                        cout << "mkdir -p " << Out_path << "/function_data_analysis_results" << endl;
                        sprintf(command, "mkdir -p %s/function_data_analysis_results", Out_path.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Make directory of function data analysis results", Error_file.c_str());

                        //Data analysis
                        cout << "sh " << path << "/Scripts/analysis/data_analysis.sh " << Out_path << "/ko_results/ko_abd.txt " << Meta_file << " " << Out_path << "/function_data_analysis_results/ dist.txt " << prefix_name << endl;
                        sprintf(command, "sh %s/Scripts/analysis/data_analysis.sh %s/ko_results/ko_abd.txt %s %s/function_data_analysis_results dist.txt %s", path.c_str(), Out_path.c_str(), Meta_file.c_str(), Out_path.c_str(), prefix_name.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Function data analysis", Error_file.c_str());
                    }
                    break;
                    
                //Multiple-species separated (Mode == 0)
                default:
                    if(species_count == 0) {
                        string error_info = "There is no suitable species for data analysis";
                        cerr << error_info << endl;
                        Echo_Error(error_info.c_str(), Error_file.c_str());
                        return 0;
                    }
                    
                    omp_set_num_threads(Min(species_count, Coren));
                    #pragma omp parallel for private(species, command)
                    for(int i = 0; i < species_count; i++) {
                        species = species_list[i];
                        
                        //Strain data analysis
                        //Mkdir
                        cout << "mkdir -p " << Out_path << "/strain_data_analysis_results/" << species << endl;
                        sprintf(command, "mkdir -p %s/strain_data_analysis_results/%s", Out_path.c_str(), species.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Make directory of strain data analysis results", Error_file.c_str());

                        //Strain data analysis for each species
                        cout << "sh " << path << "/Scripts/analysis/data_analysis.sh " << Out_path << "/strain_results/" << species << "/strain_level_abd.txt " << Meta_file << " " << Out_path << "/strain_data_analysis_results/" << species << " dist.txt " << species << endl;
                        sprintf(command, "sh %s/Scripts/analysis/data_analysis.sh %s/strain_results/%s/strain_level_abd.txt %s %s/strain_data_analysis_results/%s dist.txt %s", path.c_str(), Out_path.c_str(), species.c_str(), Meta_file.c_str(), Out_path.c_str(), species.c_str(), species.c_str());
                        outscript << command << endl;
                        Run_With_Error(command, "Strain data analysis", Error_file.c_str());

                        //Function data analysis
                        if (Is_func) {
                            //Mkdir
                            cout << "mkdir -p " << Out_path << "/function_data_analysis_results/" << species << endl;
                            sprintf(command, "mkdir -p %s/function_data_analysis_results/%s", Out_path.c_str(), species.c_str());
                            outscript << command << endl;
                            Run_With_Error(command, "Make directory of function data analysis results", Error_file.c_str());

                            //Data analysis
                            cout << "sh " << path << "/Scripts/analysis/data_analysis.sh " << Out_path << "/function_results/" << species << "/ko_abd.txt " << Meta_file << " " << Out_path << "/function_data_analysis_results/" << species << " dist.txt " << species << endl;
                            sprintf(command, "sh %s/Scripts/analysis/data_analysis.sh %s/function_results/%s/ko_abd.txt %s %s/function_data_analysis_results/%s dist.txt %s", path.c_str(), Out_path.c_str(), species.c_str(), Meta_file.c_str(), Out_path.c_str(), species.c_str(), species.c_str());
                            outscript << command << endl;
                            Run_With_Error(command, "Function data analysis", Error_file.c_str());
                        }
                    }
                    break;
            }
            //Step 4 finished

        default:
            break;
            /*
            //Parallel bdiversity R
            omp_set_num_threads(Min(command_parallel_bdiversity_scripts.size(), Coren));
            #pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < command_parallel_bdiversity_scripts.size(); i++)
            {
                char error_bdiv_parallel[BUFFER_SIZE];
                sprintf(error_bdiv_parallel, "%s/%d.log", Temp_dir.c_str(), i);
                Run_With_Error(command_parallel_bdiversity_scripts[i].c_str(), command_parallel_bdiversity_titles[i].c_str(), error_bdiv_parallel);
            }

            //Combine error
            for (int i = 0; i < command_parallel_bdiversity_scripts.size(); i++)
            {
                sprintf(command, "cat %s/%d.log >> %s", Temp_dir.c_str(), i, Error_file.c_str());
                system(command);
            }

            //Parallel R
            omp_set_num_threads(Min(command_parallel_scripts.size(), Coren));
            #pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < command_parallel_scripts.size(); i++)
            {
                char error_parallel[BUFFER_SIZE];
                sprintf(error_parallel, "%s/%d.log", Temp_dir.c_str(), i);
                Run_With_Error(command_parallel_scripts[i].c_str(), command_parallel_titles[i].c_str(), error_parallel);
            }

            //Combine error
            for (int i = 0; i < command_parallel_scripts.size(); i++)
            {
                sprintf(command, "cat %s/%d.log >> %s", Temp_dir.c_str(), i, Error_file.c_str());
                system(command);
            }
            */
    }

    if (rmdir(Temp_dir.c_str()) < 0){
        sprintf(command, "rm -rf %s", Temp_dir.c_str());
        system(command);
    }

    if (outscript) {
        outscript.close();
        outscript.clear();
    }

    Print_Report(Report_file.c_str());
    //Copy_Index(Out_path.c_str());

    cout << endl
         << "Strain2bFunc Pipeline Finished" << endl;
    outscript << endl
              << "#Strain2bFunc Pipeline Finished" << endl;
    cout << "Please check the analysis results and report at " << Out_path << endl;

    return 0;
}
