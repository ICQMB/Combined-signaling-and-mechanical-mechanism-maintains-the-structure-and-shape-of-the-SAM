//phys.h
#ifndef _PHYS_H_INCLUDED_
#define _PHYS_H_INCLUDED_
//===================
//// Forward declarations
//
////====================
//// Include dependencies
#include <math.h>
#include <random>
#include "coord.h"
////=====================
//// Simulation Constants
const double dt = .0003;
const int INIT_NUM_CYT_NODES = 15;
const int INIT_WALL_NODES = 100;
const double pi = acos(-1.0);
const double HILL_K = 1.29672531651;
const int HILL_N = 10;
const int TISSUE_BASE = 45;
//these are used to control equi angles
//and bending spring constants
const double ANGLE_FIRST_QUAD = 0.785398;
const double ANGLE_SECOND_QUAD = 2.35619;
const double ANGLE_FIRST_QUAD_Div = .43;
const double ANGLE_SECOND_QUAD_Div = 2.7;
//HIGH_ANGLE_DISCOUNT = x  means that for calculating division plane, only the bottom x% of nodes
//ordered by their angle relative to their neighbors are considered as candidates for division.
const double HIGH_ANGLE_DISCOUNT = 0.95;
//CORNER_RADIUS is the number of wall nodes away from "corners" that will.
//still be considered as "corners". 
const int CORNER_RADIUS= 4;
const double ADD_WALL_NODE_ANGLE_FIRST_QUAD =.436;
const double ADD_WALL_NODE_ANGLE_SECOND_QUAD = 2.0;
const double BOUNDARY_DAMP = 1; //TEST //was used before for boundary
const double STEM_DAMP = .1;
const double REG_DAMP = 1;
////// Cell wall mechanical parameters
const double K_BEND_STIFF = 13.5; //(13.5/15)  //29.688;//38.3622;
//17.4656;//11.8299;//52.9211;//46.4343;//4.9925;//26.4730;//35.3732;//21.2845;//6.9738;//43.5946;//48.2636;//29.6880;//39.7634;//12.7912;//31.3144;//56.4474;//58.2084;//2.1002;//20.5976;//
const double K_BEND_LOOSE =4.5433;//12.8128;
//9.9898;//2.5416;//10.9974;//11.1948;//5.2617;//6.1530;//3.8290;//1.4132;//7.5802;//6.9589;//0.8149;//4.5433;//0.6185;//9.3159;//4.8077;//2.9080;//11.8757;//8.6403;//7.9421;//
//K_LINEAR_STIFF not used
const double K_BEND_UNIFORM = 12;
const double K_LINEAR_STIFF = 0;
const double K_LINEAR_LOOSE =280.1636;//54.2730;
//157.7281;//674.3111;//511.3433;//746.0824;//461.8290;//230.2545;//134.8271;//336.2897;//581.4058;//551.7969;//84.7288;//280.1636;//688.4090;//396.0437;//612.7239;//320.9386;//381.9915;//204.3934;//477.1815;//
////Adhesion spring mechanical params
const double K_ADH = 12;
const double K_ADH_L1 = 12;
const double K_ADH_L2 = 12;
const double K_ADH_DIV = 12;
const double MembrEquLen_ADH = .9;
const double ADHThresh = 2;
const double NUMBER_ADH_CONNECTIONS = 2;
//equilibrium length of linear springs
//add node threshold length not being used
const double PERIM_INCREASE = .27;
const double MEMBR_THRESH_LENGTH = .25981; 
const double Membr_Equi_Len_Long = .07;
const double Membr_Equi_Len_Short = .07;
//
/////// Subcellular element parameters for membrane - membrane interactions
const double U_MM = 53;
const double W_MM =  0.2; //0 TEST 4
const double xsi_MM = 0.103;
const double gamma_MM = 1.1;	
//
/////// Subcellular element parameters for membrane  - internal interactions
const double U_MI = 45;
const double W_MI = 0;
const double xsi_MI = .3;
const double gamma_MI = 1.34;
const double xsi_MI_div = .4;
//
/////// Subcellular element parameters for internal - internal interactions
const double U_II = 95;
const double W_II = 6.71;
const double xsi_II = .8;
const double gamma_II = 1.34;
const double xsi_II_div = .4;


/////// Division parameters
//If the following is set to true, other division plane orientation
//rules will be overwritten to impose L1 and L2 anticlinal division.
//False will have it follow whatever rules deeper layers follow.
const bool L1_L2_FORCED_ANTICLINAL_DIV = true;

//If the following is set to true, then some of the L1 cells will "wait" for a growth phase to account for some cells growing out of the plane. 
//FREQUENTLY USED PARAMETER
//const bool OUT_OF_PLANE_GROWTH = false;

//If the following is set to true, nonlinear growth rates are applied. False
//is linear.
const bool NONLINEAR_GROWTH = true;


//If the following is set to false, then boundary division will not occur. This is also true for the L1 layer.
const bool BOUNDARY_DIVISION = false;
//All cells with layer >= STEM_LAYER will have stem damping and effects.
const int STEM_LAYER = 7;


/////// VTK parameters
//Tensile stress cytoplasm color (Calibrated from circle value)
const double CYT_COLOR = 42.0;
////=====================

//PROTOTYPE PARAMETERS

const bool OUT_OF_PLANE_GROWTH = false; //./batchGenerator -flag OOP_off
const double OOP_PROBABILITY = 0.3; //Defaults to 0.3
const int Weird_WUS = 0;
const int Div_ON;


//extern int NUM_STEPS_PER_FRAME
//extern int RECENT_DIV_NUM_FRAMES
#endif

