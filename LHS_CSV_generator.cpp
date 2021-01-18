#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <random>
//
using namespace std;
double unifRand();
int unifRandInt(int a, int b);
double unifRand(double a, double b);
double stratifiedBins(double a, double b, int k, int N);
void swapEntries(vector<int>& X, unsigned int a, unsigned int b);
void seed();

//Produces LHS vectors into a CSV for performing
//Sobol' GSA. 

int main(int argc, char* argv[]) { 

	//seed random generator
	mt19937::result_type seed = time(0);

	//Get output name
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
	//Parameter Data Type
	//1: double
	//2: int
	//3: bool
	vector<int> param_types;
	param_types.clear();
	cout << endl << "Enter data types of the parameters: "
		<< endl
		<< "1: double" << endl << "2: int" << endl << "3: boolean" << endl;
	int temp_int;
	for (unsigned int i = 0; i < P; i++) {
		cout << "Parameter " << i+1 << ": ";
		cin >> temp_int;
		if (temp_int < 1 || temp_int > 3) { 
			cout << "Invalid entry. Exiting..." << endl;
			exit(1);
		}
		param_types.push_back(temp_int);
	}
	vector<string> minvals;
	vector<string> maxvals;
	//Get Parameter Ranges
	for (unsigned int i = 0; i < P; i++) { 
		if (param_types.at(i) == 1 || param_types.at(i) == 2) { 
			cout << "Enter parameter " << i+1 << " min:";
			cin >> temp_str;
			minvals.push_back(temp_str);
			cout << "Enter parameter " << i+1 << " max:";
			cin >> temp_str;
			maxvals.push_back(temp_str);
		} else { 
			minvals.push_back("0");
			maxvals.push_back("1");
		}

	}
	for (unsigned int i = 0; i < P; i++) { 
		cout << "Parameter " << i + 1 << " info:" << endl;
		cout << "Minval: " << minvals.at(i) << endl;
		cout << "Maxval: " << maxvals.at(i) << endl;
		cout << "Type: " << param_types.at(i) << endl;
	}
	
	//Get number of samples
	int num_samples;
	cout << "Enter number of samples:" << endl;
	cin >> num_samples;
	if (num_samples < P) { 
		cout << "number of samples must be >= P" << endl;
		exit(1);
	}
	vector<vector<int>> X_sample_bins;
	vector<vector<int>> Y_sample_bins;

	for (int i = 0; i < P; i++) { 
		X_sample_bins.push_back(vector<int>());
		Y_sample_bins.push_back(vector<int>());
	}
	for (int k = 0; k < P; k++) {
		for (int i = 0; i < num_samples; i++) { 
			X_sample_bins.at(k).push_back(i);
			Y_sample_bins.at(k).push_back(i);
		}
		for (int j = 0; j < num_samples; j++) { 
			swapEntries(X_sample_bins.at(k),j,unifRandInt(j,num_samples-1));
			swapEntries(Y_sample_bins.at(k),j,unifRandInt(j,num_samples-1));
		}
	}
	double start_bracket;
	double end_bracket;
	double bool_roll;
	int int_roll;
	double double_roll;

	vector<string> X_k;
	vector<string> Y_k;
	
	char delimiter;
	
	for (unsigned int i = 0; i < num_samples; i++) { 
		//Choose bracket to sample; try again if already sampled
	 	//Generate X and Y 
		X_k.clear();
		Y_k.clear();
		for (unsigned int j = 0; j < P; j++) { 
			switch(param_types.at(j)) { 
				case 1: 
					//Makes uniform rand double into X and Y 
					//X
					start_bracket = stod(minvals.at(j));
					end_bracket = stod(maxvals.at(j));
					double_roll = stratifiedBins(
							start_bracket,
							end_bracket,
							X_sample_bins[j][i],
							num_samples);
					X_k.push_back(to_string(double_roll));

					//Y
					double_roll = stratifiedBins(
							start_bracket,
							end_bracket,
							Y_sample_bins[j][i],
							num_samples);
					Y_k.push_back(to_string(double_roll));
										
					break;
				case 2: //int
					//Make uniform rand int in range.
					//Makes uniform rand int into X and Y 
					//X
					start_bracket = stoi(minvals.at(j));
					end_bracket = stoi(maxvals.at(j)) + 1;
					int_roll = 
						floor ( 
						stratifiedBins(
							start_bracket,
							end_bracket,
							X_sample_bins[j][i],
							num_samples
							)
						);
					if(int_roll < stoi(minvals.at(j))){
						int_roll++;
					} else if (int_roll > stoi(maxvals.at(j))) { 
						int_roll--;
					}
					X_k.push_back(to_string(int_roll));


					//Y

					int_roll = 
						floor ( 
						stratifiedBins(
							start_bracket,
							end_bracket,
							Y_sample_bins[j][i],
							num_samples
							)
						);
					if(int_roll < stoi(minvals.at(j))){
						int_roll++;
					} else if (int_roll > stoi(maxvals.at(j))) { 
						int_roll--;
					}
					Y_k.push_back(to_string(int_roll));
					break;
				case 3: //bool
					start_bracket = 0;
					end_bracket = 2;
					bool_roll = stratifiedBins(
								start_bracket,
								end_bracket,
								X_sample_bins[j][i],
								num_samples);
					if (bool_roll <= 1) { 
						X_k.push_back("0");
					} else { 
						X_k.push_back("1");
					}
					//Y
					bool_roll = stratifiedBins(
								start_bracket,
								end_bracket,
								Y_sample_bins[j][i],
								num_samples);
					if (bool_roll <= 1) { 
						Y_k.push_back("0");
					} else { 
						Y_k.push_back("1");
					}
					break;
				default:
					cout << "Unrecognized data type";
					exit(1);
			}
		} 


		//Generate CSV that has the structure
		//X
		//Y
		//X_1 (= Y_-1)
		//...
		//X_z (= Y_-z)
		//...
		//X_P (= Y_-P)

		//Print X to CSV
		for (int j = 0; j < P; j++) { 
			//Delimiter line makes the last cell end with newline 
			//instead of comma.
			delimiter = (j == P-1) ? '\n' : ','; 
			CSV << X_k.at(j) << delimiter;
		}
		//Print Y to CSV
		for (int j = 0; j < P; j++) { 
			delimiter = (j == P-1) ? '\n' : ','; 
			CSV << Y_k.at(j) << delimiter;
		}
		//Print X_z, or Y_-z to CSV for all parameter values.
		for (int z = 0; z < P; z++) { 
			for (int j = 0; j < P; j++) { 
				delimiter = (j == P-1) ? '\n' : ','; 
				if (j == z) {
					CSV << Y_k.at(j) << delimiter;
				} else { 
					CSV << X_k.at(j) << delimiter;
				}
			}
		}
	}
	return 0;
}

//Generates N bins out of continuous [a,b] and selects
//randomly from the kth out of N equal bins. Bins
//are labeled bin 0,1,..., N-2, N-1.
double stratifiedBins(double a, double b, int k, int N) { 
	
	if (b <= a) { 
		cout << "Problem! b < a in requested bins." << endl;
		exit(1);
	}
	double lower_bracket = a + ((b-a)/static_cast<double>(N))*static_cast<double>(k);
	double upper_bracket = a + ((b-a)/static_cast<double>(N))*static_cast<double>(k+1);
	return unifRand(lower_bracket,upper_bracket);
}

//swaps elements a and b in integer vector x.
void swapEntries(vector<int>& X, unsigned int a, unsigned int b) {
	if (a >= X.size() || b >= X.size()) { 
		cout << "Trying to access outside of X!" << endl;
		exit(1);
	} else if (a == b) return;
	int temp;
	temp = X.at(a);
	X.at(a) = X.at(b);
	X.at(b) = temp;
	return;
}

int unifRandInt(int a,int b) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(a,b);
	return dis(gen);	
}
double unifRand() {
	return rand() / double (RAND_MAX);
}
double unifRand(double a, double b) {
	return (b-a)*unifRand() + a;
}
