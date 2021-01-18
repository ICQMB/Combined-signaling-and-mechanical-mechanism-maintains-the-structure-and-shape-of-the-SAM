//
//cell_calibration.cpp
//========================
//Forward Declarations
//


//========================
//Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//==========================

//=====================================================
//============================
// Getters and Setters 
// ===========================
// ====================================================

//this function returns vector of strain measurements
//stored at each time point during stretching
void Cell::get_Strain(vector<double>& strain) {
	strain = this->strain_vec;
	return;
}

//this function returns vector of stress measurements
//stored at each time point during stretching
void Cell::get_Stress(vector<double>& stress) {
	stress = this->stress_vec;
	return;
}

//======================================================
//==================================
//Functions for calibration
//==================================
//======================================================
//measuring pressure
double Cell::compute_pressure(){
	Wall_Node* curr = this->top;
	Wall_Node* next = NULL;
	Coord force = Coord(0,0);
	double curr_length = 0;
	double total_length = 0;
	double pressure = 0;
	force += top->get_CytForce();

	for(unsigned int i=0; i < counter_left-1; i++) {
		next = curr->get_Left_Neighbor();
		force += next->get_CytForce();
		curr_length = (curr->get_Location() - next->get_Location()).length();
		total_length += curr_length;
		curr = next;
	}
	for(unsigned int i=0; i < counter_right-1; i++) {
		next = curr->get_Right_Neighbor();
		force += next->get_CytForce();
		curr_length = (curr->get_Location() - next->get_Location()).length();
		total_length += curr_length;
		curr = next;
	}
	
//	cout << "Area: " << this->calc_Area() << endl;
	pressure = force.length()/total_length;
	
	return pressure;
	cout << "Pressure" << pressure << endl;
}

//compression test functions
void Cell::set_Stationary_Points(int Ti) {
	if(Ti == calibStart) {
		this->top = this->closest_node_top();
		this->top->set_isStationary();
		this->bottom = this->closest_node_bottom();
		this->bottom->set_isStationary();
	}
	Wall_Node* curr = this->top;
	Wall_Node* next = NULL;
	double y_coord = top->get_Location().get_Y();
	do {
		next = curr->get_Left_Neighbor();
		curr = next;
	} while (next->get_isStationary());

	if((abs(next->get_Location().get_Y() - y_coord) < .02)){
			next->set_isStationary();
	}
	
	curr = this->top;
	do {
		next = curr->get_Right_Neighbor();
		curr = next;
	} while (next->get_isStationary());
	
	if(abs(next->get_Location().get_Y() - y_coord) < .02){
			curr->set_isStationary();
	}
	curr = this->bottom;
	next = NULL;
	y_coord = bottom->get_Location().get_Y();
	do {
		next = curr->get_Left_Neighbor();
		curr = next;
	} while (next->get_isStationary());

	if((abs(next->get_Location().get_Y() - y_coord) < .02)){
			next->set_isStationary();
	}
	
	curr = this->bottom;
	do {
		next = curr->get_Right_Neighbor();
		curr = next;
	} while (next->get_isStationary());
	
	if(abs(next->get_Location().get_Y() - y_coord) < .02){
			curr->set_isStationary();
	}
			
	return;
}

//elastic modulus/stretching test functions
void Cell::compress(){
	//stretch the cell
	Wall_Node* curr = this->top;
	Wall_Node* next = NULL;
	this-> counter_left = 0;
		do {
			next = curr->get_Left_Neighbor();
			curr->pull_node();
			counter_left++;
//			cout << "pulling" << endl;
			curr = next;
			
		} while (next->get_isStationary());

	curr = top;
	this->counter_right = 0;
		do {
			next = curr->get_Right_Neighbor();
			curr->pull_node();
			counter_right++;
//			cout << "pulling" << endl;
			curr = next;
		} while(next->get_isStationary());
	return;
}
/*void Cell::extensional_strain() {
    double new_area = this->calc_Area();
	double delta_L = new_area -this-> curr_area;
	double strain = delta_L/curr_area;
	this->curr_area = new_area;
	strain_vec.push_back(strain);
	return;
}

void Cell::tensile_Stress(){
	Wall_Node* curr = this->left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	Coord force_sum = Coord(0,0);
	double force;
	double curr_length = 0;
	double total_length = 0;

	do {
		force_sum += curr->get_f_EXT();
		next = curr->get_Left_Neighbor();
		curr_length = (curr->get_Location() - next->get_Location()).length();
		total_length += curr_length;
		curr = next;
	} while (next != orig);

	force = force_sum.length();
	stress_vec.push_back(force/total_length);
	return;
 }*/


