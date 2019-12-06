#ifndef GRID_HPP
#define GRID_HPP

#include "./LSH_Structure.hpp"

using namespace std;

class Point;
class Curve;
class NN;
class LSH;

class Grid{
		vector<double> a;	// beginning of the grid
		pair<double, double> convert_point(pair<double, double>);
	public:
		Grid(double delta);
		vector<double> snap(Curve* c, double* max = NULL);
};

class Grid_LSH{
		int max_d;
		vector<Curve*> curveset;
		vector<Grid*> Grids;
		vector<double> snap_num;
		vector<LSH*> lsh;
	public:
		Grid_LSH(vector<Curve*>* curves, int L_grid, int k_vec, int max_d, int min_d);
		~Grid_LSH();
		NN* predict(Curve* curve, int max_d);
		vector <Curve *> * get_bucket(Curve * curve);
};

// class Grid_Hypercube{

// 		vector<Curve*> curveset;
// 		vector<Grid*> Grids;
// 		vector<double> snap_num;
// 		vector<Hypercube*> hypercube;
// 	public:
// 		Grid_Hypercube(vector<Curve*>* curves, unsigned int k_hypercube, unsigned int M, unsigned int probes, int L_grid, int max_d, int min_d);
// 		~Grid_Hypercube();
// 		NN* predict(Curve* curve, int max_d);
// };

#endif
