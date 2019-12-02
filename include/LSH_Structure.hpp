#ifndef LSH_STRUCTURE_HPP
#define LSH_STRUCTURE_HPP

#include <random>
#include <unordered_map>
#include <vector>
#include "./utilities.hpp"

using namespace std;

class Point;
class Curve;
class NN;

class H{
		vector<double> s;
		int m;
		double w;
		int M;

		int modulo(int a, int b);
		int modPow(int b, int e, int m);
	public:
		H(int dimension, int m, double w, int M);
		inline ~H() { this->s.clear(); }
		int operator[](Point* point);	// Returns the h value of the point
};

class G{
		vector<H*> h;
	public:
		G(int k, int dimension, int m, double w);
		~G();
		unsigned int operator[](Point* point);	// Returns the g value of the point
};

class LSH{
		vector<Point*> pointset;
		vector<G*> g;
		vector<unordered_map<unsigned int, vector<Point*>>> hash_tables;

		void train();
	public:
		LSH(vector<Point*>* points, int L, int k, int dimension);
		~LSH();
		NN* predict(Point * point, int r = 0);
};

#endif