//Include Dependencies
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//REQUIRES UNIX ENGINE
#include <unistd.h>
#include <bits/stdc++.h>

using namespace std;

//File accepts a vector of parameters, format:
//First line is the simulation flags followed by their values.
//
//Second line is formatted as flag_1, ... ,flag_p
//where flags are those included for main to run the file.
//Third line onwards are the actual vectors of parameters to be fed into 
//main and batch.
int main(int argc, char* argv[]) { 

	//Get file or exit.
	if (argc == 1) { 
		cout << "Please provide CSV file name" << endl; 
		return 0;
	}
	string CSVname = argv[1];
	//Declare ifstream CSV
	ifstream CSV(CSVname.c_str()); 
	if(!CSV) { 
		cout << CSVname.c_str()  << " is not available" << endl;
		return 0;
	}


	//initialize variables
	string cell;
	string line;	
	//Iterator for the loop
	int k;
	cout << "Please enter integer k for row 1's index: ";
	cin >> k;
	k = k-1;
	vector<string> batch_info;
	batch_info.clear();
	vector<string> parameter_names;
	parameter_names.clear();
	vector<string> parameters;
	parameters.clear();
	stringstream command;


	//Change this to be the name of the .out file for generating batch script.
	string command_base = "./batchGenerator.out ";
	getline(CSV,line);
	stringstream ss(line);
	//Get parameter names to store flags, and batch info.
	//test_name_index
	unsigned int NI = 0;
	unsigned int test_name_index;
	while(getline(ss,cell,',')) { 
		NI++;
		if (!strcmp(cell.c_str(),"-test")) { 
			test_name_index = NI;
		}
		batch_info.push_back(cell);
	}
	if (!ss && cell.empty()) {
		batch_info.push_back("");
	} 
	getline(CSV,line);
	ss.str("");
	ss.clear();
	ss.str(line);
	while(getline(ss,cell,',') ){ 
		parameter_names.push_back(cell);	
	}
	//Loop through the kth row, build and submit the simulation.
	while(getline(CSV,line)) { 
		k++;
		parameters.clear();
		ss.str("");
		ss.clear();
		ss.str(line);
		while(getline(ss,cell,',')) {
			parameters.push_back(cell);
		}
		cout << "Parameter set "  << k << endl;
		for (unsigned int i = 0; i < parameters.size(); i++) { 
			cout << parameters.at(i) << " " << flush;
		}
		//Generate batch file with parameternames, parameters
		//Clears the "command" string stream.
		command.str("");
		command.clear();
		command << command_base;
		for (unsigned int i = 0; i < batch_info.size(); i++) {
			if (i == test_name_index) { 
				command << batch_info.at(i) << "_" << k << " ";
			} else { 
				command << batch_info.at(i) << " ";
			}
		}
		for (unsigned int i = 0; i < parameters.size(); i++) { 
			command << "-par " + parameter_names.at(i) + " " + parameters.at(i) + " ";
		}
		cout << endl << "Command: " << command.str() << endl;
		//Submits command to bash to generate AUTO_BATCH.sh.
		system(command.str().c_str());
		//Wait 2 seconds for batch file to generate.	
		sleep(2);
		//Submit batch file
		system("sbatch AUTO_BATCH.sh");
		//Wait 6 seconds for batch to submit properly.
		sleep(6);
	}
	return 0;
}
