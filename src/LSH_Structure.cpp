#include "../include/LSH_Structure.hpp"

using namespace std;

H::H(int dimension, int m, double w, int M){

	// Create s
	random_device rd;
	mt19937 generator(rd());
	uniform_real_distribution<> distribution(0.0,w/1.0);

	for (int i = 0; i < dimension; ++i)
	{
		this->s.push_back(distribution(generator));
	}

	this->m = m;
	this->w = w;
	this->M = M;
}

// Returns the h value of the point
int H::operator[](Point* point){

	int sum = 0;
	int a = 0;
	int d = point->get_dimension();
	for (int i = 0; i < d; i++)
	{
		a = modulo((floor(((*point)[i] - this->s[i])/(this->w)))*(modPow(this->m, d-i-1, this->M)), this->M);
		sum += a;
	}
	return modulo(sum, this->M);
}

int H::modulo(int a, int b) {
  int m = a % b;
  if (m < 0) {
    m = (b < 0) ? m - b : m + b;
  }
  return m;
}

int H::modPow(int b, int e, int m){

	int r = 1;
	for (int i = 0; i < e; ++i)
	{
		r = (r*b)%m;
	}
	return r;
}

G::G(int k,int dimension, int m, double w){

	//Create h
	for (int i = 0; i < k; ++i)
	{
		try{
			this->h.push_back(new H(dimension, m, w, pow(2, floor(32/k))));
		}
		catch(std::bad_alloc&){

			for (int j = 0; j < i; ++j)
			{
				delete this->h[j];
			}
		}
	}
}

G::~G(){

	delete_vector<H>(&(this->h));
	this->h.clear();
}

// Returns the g value of the point
unsigned int G::operator[](Point* point){

	unsigned int g = 0;
	int k = this->h.size();
	int shift = floor(32/k);

	// Concatenate h values
	for (int i = 0; i < k; ++i)
	{
		g = g << shift;
		g = g | (*(this->h[i]))[point];
	}
	return g;
}


LSH::LSH(vector<Point*>* points, int L, int k, int dimension){

	// Create pointset
	for (unsigned int i = 0; i < (*points).size(); ++i)
	{
		this->pointset.push_back((*points)[i]);
	}

	double w = 4*average_distance(&(this->pointset)) + 7;
	int m = 3;

	// Create the g functions and the hash tables
	for (int i = 0; i < L; ++i)
	{
		this->g.push_back(new G(k, dimension, m, w));
		unordered_map<unsigned int, vector<Point*>> hash_table;
		this->hash_tables.push_back(hash_table);
	}

	// Fill the hashtables
	this->train();
}

LSH::~LSH(){

	delete_vector<Point>(&(this->pointset));
	this->pointset.clear();
	delete_vector<G>(&(this->g));
	this->g.clear();
	for (unsigned int i = 0; i < this->hash_tables.size(); ++i)
	{
		unordered_map<unsigned int, vector<Point*>>:: iterator itr; 
	    for (itr = this->hash_tables[i].begin(); itr != this->hash_tables[i].end(); itr++) 
    	{ 
	        itr->second.clear(); 
    	}
	}
	this->hash_tables.clear();
}

void LSH::train(){

	int L = (this->g).size();

	// Fill the hash tables with the pointset
	unsigned int key = 0;
	for (unsigned int j = 0; j < this->pointset.size(); ++j)	// For every point
	{
		for (int i = 0; i < L; ++i)	// For every g
		{
			key = (*((this->g)[i]))[(this->pointset)[j]];
			if ( this->hash_tables[i].find(key) == this->hash_tables[i].end() )	// If the key does not exist
			{
				// Insert the key
				vector<Point*> p;
				p.push_back(this->pointset[j]);
				this->hash_tables[i].insert(make_pair(key, p));
			}
			else	// If the key exists
				this->hash_tables[i].at(key).push_back(this->pointset[j]);
		}		
	}

}

// Find nearest neighbor
NN* LSH::predict(Point* point, int r){

	int min_distance, distance;
	string min_id = "";
	set<string> near_neighbors;

	int L  = (this->g).size();

	for (int i = 0; i < L; ++i) // For every g
	{
		unsigned int key = (*(this->g[i]))[point];	// Find g value of the point

		if ( (this->hash_tables)[i].find(key) != (this->hash_tables)[i].end() ) // If the key exists
		{
			vector<Point*> buckets_points = (this->hash_tables)[i].at(key);
			Point* p;

			int size = buckets_points.size();
			
			int step = 1;

			if ( r != 0 ){
				// Increase step according to bucket size
				step = size/500 + 1;
			}

			for (int j = 0; j < size; j+=step)	// For every point in the bucket
			{	
				p = buckets_points[j];
				distance = manhattan_dist(p, point);	// Calc the distance between the 2 points
				
				// Update near neighbors
				if ( distance < r )
					near_neighbors.insert(p->get_id());

				// Update nearest neighbor
				if ( min_id == "" )	// If this is the first point
				{
					min_id = p->get_id();
					min_distance = distance;
				}
				else
				{
					if ( distance < min_distance )
					{
						min_distance = distance;
						min_id = p->get_id();
					}
				}
			}
		}
	}

	// If i have found a nearest neighbor
	if ( min_id.compare("") !=0 )
	{
		near_neighbors.erase(min_id);

		NN * nearest_neighbor;
		try{
			nearest_neighbor = new NN(min_id, min_distance, &near_neighbors);
		}
		catch(bad_alloc&)
		{
			cerr << "No memory available" << endl;
			return NULL;
		}

		return nearest_neighbor;
	}
	else
		return NULL;
}
