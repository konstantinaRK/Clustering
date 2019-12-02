#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <bits/stdc++.h> 
#include "./classification.hpp"

using namespace std;

class Classification;

class Point{
		string id;
		vector<double> X;
	public:
		Point(string id, vector<double>* X);
		~Point();
		inline string get_id(void) { return this->id; }
		int get_dimension(void);
		double operator[](unsigned int i);
};

class Curve
{
	private:
		string curve_id;
		vector <pair <double, double>> points;
	public:
		Curve(string id);
		~Curve();
		inline string get_id(void){ return this->curve_id; }
		void add_point(double x, double y);
		int get_length(void);
		pair<double, double> operator[](int i);
		void clear(void);
	
};

class NN {
		string id;
		double distance;
		vector<string> r_near_neighbors;
	public:
		NN(string id, double distance, set<string>* neighbors=NULL);
		~NN();
		inline string get_id(void){ return this->id; } 
		inline double get_distance(void){ return this->distance; }
		unsigned int r_near_neighbors_size();
		string get_near_neighbor(int i);
		vector<string> get_near_neighbors(void);
		void add_neighbors(vector<string>* new_neighbors);
};


Classification * dataHandling(int argc, char * argv[], string * output_file, bool * complete);


template <class C>
void delete_vector(vector<C*>* v){

	for (unsigned int i = 0; i < (*v).size(); ++i)
	{
		delete (*v)[i];
	}

	(*v).clear();
}

double min(double x, double y, double z, string* direction);

double DTW_distance(Curve* x1, Curve* x2, vector<pair<int, int>>* opt_trav=NULL);
double eucl_dist(pair<double, double> x, pair<double, double> y);
double manhattan_dist(Point* x, Point* y);
double average_distance(vector<Point*>* pointset);
NN* brute_force(Point* point, vector<Point*>* pointset);

bool read(string file_name, vector<Point*>* points);
bool point_proccessing(vector<Point*>* points, string p, int d = -1);


// TODO
double DTW(Curve *x1, Curve * x2);

#endif