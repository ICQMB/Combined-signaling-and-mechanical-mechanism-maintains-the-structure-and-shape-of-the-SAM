#include <iostream>
#include <cstring>
#include <fstream>
#include <math.h>
#define PI 3.14159265
using namespace std;

//Makes a file "18rings.txt" that represents an initial condition of 18 concentric circles

int main(int argc, char* argv[]) {

	ofstream ofs;
	ofs.open("staggered_generated.txt");
	int Layers = 8;
	int N = 0;
	const double INIT_RADIUS = 3.75;
	double lamella_width_horiz = 0.3;
	double lamella_width_vert = 0.3;
	double offset;
	double center_spacing_h, center_spacing_v;
	double x, y;
	double layer_start;
	for (int i = 1; i < argc; i++) { 
		if (!strcmp(argv[i], "-h")) { 
			lamella_width_horiz = stod(argv[i+1]);
		} else if (!strcmp(argv[i], "-v")) { 
			lamella_width_vert = stod(argv[i+1]);
		}
	}
	center_spacing_h = INIT_RADIUS * 2 + lamella_width_horiz;
	center_spacing_v = INIT_RADIUS * 2 + lamella_width_vert;

	for (int i = 1; i <= Layers; i++) {
		int cells_this_layer = Layers+4-i;
		if (i == Layers) {
			cells_this_layer++;
		}
		y = (1-i) * center_spacing_v;
		offset = center_spacing_h / 2.0;
		if (cells_this_layer%2 == 1) { 
			ofs << "CellRank:" << N << '\n';
			ofs << "Center:" << 0 << "," << y << '\n';
			ofs << "Radius:" << INIT_RADIUS << "\n";
			ofs << "Layer:" << i << "\n";
			ofs << "Boundary:0\n";
			ofs << "End_Cell\n\n";
			N++;
			cells_this_layer--;
			offset = 0;

		}
		for (int j = 1; j <= cells_this_layer/2; j++)
		{
			x = (j) * center_spacing_h - offset;
			ofs << "CellRank:" << N << '\n';
			ofs << "Center:" << x << "," << y << '\n';
			ofs << "Radius:" << INIT_RADIUS << "\n";
			ofs << "Layer:" << i << "\n";
			ofs << "Boundary:";
			if (j == cells_this_layer / 2) { 
				ofs << 1;
			} else {
				ofs << 0;
			}
			ofs << "\nEnd_Cell\n\n";
			N++;

			ofs << "CellRank:" << N << '\n';
			ofs << "Center:" << (-1) * x << "," << y << '\n';
			ofs << "Radius:" << INIT_RADIUS << "\n";
			ofs << "Layer:" << i << "\n";
			ofs << "Boundary:";
			if (j == cells_this_layer / 2) { 
				ofs << 1;
			} else {
				ofs << 0;
			}
			ofs << "\nEnd_Cell\n\n";
			N++;
		}
	}
	ofs.close();


	return 0;
}
