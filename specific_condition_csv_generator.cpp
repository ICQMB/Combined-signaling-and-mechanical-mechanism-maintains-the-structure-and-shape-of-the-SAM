#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <random>
//
using namespace std;

int main(int argc, char* argv[]) { 

	//Get output name
	int num_reps;
	string CSVname;
	cout << "Please enter output CSV file name (no file extension): " << endl;
	cout << ">";
	cin >> CSVname;
	
	//Generate CSV file
	ofstream CSV((CSVname + ".csv").c_str());

	//Get Batch Info from argv[] 
	cout << "Enter batch info: Comma separated, no trailing comma, no spaces:" 
		<< endl << ">";
	string temp_str;
	cin >> temp_str;
	CSV << temp_str << endl;
	//Get number of Parameters / Parameter Flags
	cout << "Enter parameter flags: Comma separated, no trailing commas, no spaces:" << endl << ">";
	cin >> temp_str;
	CSV << temp_str << endl;
	unsigned int P;
	cout << "Number of parameters being varied?" << endl;
	cout << ">";
	cin >> P;
	int temp_int;	
	cout << "Number of repetitions?" << endl;
	cout << ">";
	cin >> num_reps;
	cout << endl << "How many parameter values per variable?: ";
	vector<int> num_values_per_param;
	for (unsigned int i = 0; i < P; i++) {
		cout << "Parameter " << i+1 << ": ";
		cin >> temp_int;
		num_values_per_param.push_back(temp_int);
	}
	vector<vector<string>> paramvals;
	
	//Get Parameter Ranges
	for (unsigned int i = 0; i < P; i++) { 
		vector<string> param_i_vals;
		for (int j = 0; j < num_values_per_param.at(i); j++) { 
			cout << "Enter parameter " << i+1 << " value " << j+1 << ":";
			cin >> temp_str;
			param_i_vals.push_back(temp_str);
		}
		paramvals.push_back(param_i_vals);
	}
	
	//Get number of samples
	int num_samples = 1;
	for (unsigned int i = 0; i < num_values_per_param.size(); i++) { 
		num_samples *= num_values_per_param.at(i);
	}

	vector<string> X;

	char delimiter;

	for (unsigned int i = 0; i < num_samples; i++) { 
		//Choose bracket to sample; try again if already sampled
		//Generate X and Y 
		vector<int> entries_to_access;
		for (unsigned int k = 0; k < P; k++) { 
			entries_to_access.push_back(0);
		}
		for (int j = 1; j <= i; j++) {
			entries_to_access.at(0) = entries_to_access.at(0) + 1;
			for (int k = 0; k < P; k++) { 
				if (entries_to_access.at(k) >= num_values_per_param.at(k)) {
					entries_to_access.at(k) = 0;
					entries_to_access.at(k+1) = entries_to_access.at(k+1) + 1;
				}
			}
		}
		for (int r = 0; r < num_reps; r++) { 
			for (int j = 0; j < P; j++) { 
				//Delimiter line makes the last cell end with newline 
				//instead of comma.
				delimiter = (j == P-1) ? '\n' : ','; 
				CSV << paramvals.at(j).at(entries_to_access.at(j)) << delimiter;
			}
		}
	}
	return 0;
}
