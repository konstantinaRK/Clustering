#ifndef CLASSIFICATION_HPP
#define CLASSIFICATION_HPP

#include <vector>
#include <set>
#include "./utilities.hpp"

using namespace std;

class Point;
class Curve;
class Clustering;

class Classification
{
	protected:
		int cluster_num;
		int grid_num;
		int vector_htables_num;
		int vector_hfunc_num;
		vector<Clustering *> clusterings;
	public:
		Classification(string config);
		virtual ~Classification();
		// void train();
			
};

class Classification_Points: public Classification
{
	private:
		vector <Point *> data;
	public:
		Classification_Points(string input_file, string config, short int flag);
		~Classification_Points();
	
};

class Classification_Curves: public Classification
{
	private:
		vector <Curve *> data;
	public:
		Classification_Curves(string input_file, string config, short int flag);
		~Classification_Curves();
	
};

class Clustering
{
	protected:
		set <int> centers;
		short int flag;
		multimap<int, int> clusters;
	public:
		// Clustering(short int flag, int cluster_num, vector<Point*>* data){};
		virtual ~Clustering(){};
		virtual void initialization1(int cluster_num, int data_size);
		virtual void initialization2(){};
		virtual void assignment1(){};
		virtual void assignment2(){};
		virtual void update1(){};
		virtual void update2(){};
	
};

class Point_Clustering: public Clustering
{
	private:
		double binary_search(vector<double>* P, double x);
		double min_dist(vector<Point*>* data, int pos);
	public:
		Point_Clustering(short int flag, int cluster_num, vector<Point*>* data);
		// ~Point_Clustering(){};
		void initialization2(int cluster_num, vector<Point*>* data);
		void assignment1(vector<Point*>* data);
		void assignment2(){};
		void update1(){};
		void update2(){};
	
};

class Curve_Clustering: public Clustering
{
	private:
	public:
		Curve_Clustering(short int flag, int cluster_num, vector<Curve*>* data);
		// ~Curve_Clustering(){};
		void initialization2(){};
		void assignment1(){};
		void assignment2(){};
		void update1(){};
		void update2(){};
	
};

#endif