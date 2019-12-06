#ifndef CLASSIFICATION_HPP
#define CLASSIFICATION_HPP

#include <vector>
#include <set>
#include <utility>
#include <unordered_map>
#include "./utilities.hpp"
#include "./Grid.hpp"
#include "./LSH_Structure.hpp"

using namespace std;

class Point;
class Curve;
class Clustering;
class Grid_LSH;
class LSH;

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
};

class Classification_Points: public Classification
{
	private:
		vector <Point *> data;

		LSH * lsh;

	public:
		Classification_Points(string input_file, string output_file, string config, short int flag, bool complete);
		~Classification_Points();
	
};

class Classification_Curves: public Classification
{
	private:
		vector <Curve *> data;

		Grid_LSH * grid_lsh;

	public:
		Classification_Curves(string input_file, string output_file, string config, short int flag, bool complete);
		~Classification_Curves();
	
};

class Clustering
{
	protected:
		vector <bool> mean_centers;
		short int flag;
		multimap<int, int> clusters;

		template<class D>
		double min_dist(vector<D*>* data, vector<D*>* centers, int pos);

		template<class D>
		bool centers_changed(vector<D*>* old_centers, vector<D*>* cur_vectors, vector<D*>* data);

		template<class D>
		int calc_2min_cen(D* d, vector<D*>* centers);

		template<class D>
		double cluster_dist(D* d, int center_pos, vector<D*>* data);

		template<class D>
		double Cluster_Silhouette(vector<D*>* data, vector<D*>* centers, int cluster_num);
	public:
		Clustering(int cluster_num);
		virtual ~Clustering();

		template<typename vector_type>
		void initialization1(unsigned int cluster_num, vector <vector_type*> * centers, vector<vector_type*>* data);
		template<class D>
		void initialization2(int cluster_num, vector<D*>* centers, vector<D*>* data);
		template<class D>
		void assignment1(vector<D*>* centers, vector<D*>* data);
		template <typename lsh_type, typename vector_type>
		void assignment2(lsh_type * lsh, vector <vector_type *> * data, vector <vector_type *> *centers);
		template<class D>
		bool update1(vector<D*>* centers, vector<D*>* data);

		virtual bool update2(vector<Point*>* data){return false;};

		virtual double distance(Point *, Point *){return 0;};
		virtual double distance(Curve *, Curve *){return 0;};

		virtual void write_output(string output_file, double time, bool optional, bool means, vector<Point*>* data){};
		virtual void write_output(string output_file, double time, bool optional, bool means, vector<Curve*>* data){};

};

class Point_Clustering: public Clustering
{
	private:
		vector <Point *> centers;
	public:
		Point_Clustering(short int flag, int cluster_num, vector<Point*>* data, LSH* lsh);
		~Point_Clustering();
		
		bool update2(vector<Point*>* data);

		double distance(Point *c1, Point *c2);
		void write_output(string output_file, double time, bool optional, bool means, vector<Point*>* data);
};

class Curve_Clustering: public Clustering
{
	private:
		vector <Curve *> centers;
	public:
		Curve_Clustering(short int flag, int cluster_num, vector<Curve*>* data, Grid_LSH* grid_lsh, int min_d, int max_d);
		~Curve_Clustering();
		
		bool update2(vector<Curve*>* data);
	
		double distance(Curve *c1, Curve *c2);
		void write_output(string output_file, double time, bool optional, bool means, vector<Curve*>* data);
};

#endif