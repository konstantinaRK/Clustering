#include "../include/classification.hpp"

#include <fstream>
#include <random>
#include <time.h>

// Class Classification functions
Classification::Classification(string config)
{
	#if DEBUG
		cout << "Classification in" << endl;
	#endif

	this->cluster_num = -1;
	this->grid_num = 2;
	this->vector_htables_num = 2;
	this->vector_hfunc_num = -1;

	ifstream vars;
	vars.open(config);
	if (!vars.is_open())
	{
		throw invalid_argument("Wrong config file name.");
	}

	// Read values and initialize variables
	string line;
	while (getline(vars, line))
	{
		if (!line.compare(0, 19, "number_of_clusters:"))
		{
			line = line.substr(line.find(":")+1);
			this->cluster_num = stoi(line);
			if (this->cluster_num < 1)
			{
				throw invalid_argument("Invalid cluster number value.");
			}
		}
		else if (!line.compare(0, 16, "number_of_grids:"))
		{
			line = line.substr(line.find(":")+1);
			this->grid_num = stoi(line);
			if (this->grid_num < 1)
			{
				throw invalid_argument("Invalid grid number value.");
			}
		}
		else if (!line.compare(0, 29, "number_of_vector_hash_tables:"))
		{
			line = line.substr(line.find(":")+1);
			this->vector_htables_num = stoi(line);
			if (this->vector_htables_num < 1)
			{
				throw invalid_argument("Invalid vector hash tables value.");
			}
		}
		else if (!line.compare(0, 32, "number_of_vector_hash_functions:"))
		{
			line = line.substr(line.find(":")+1);
			this->vector_hfunc_num = stoi(line);
			if (this->vector_hfunc_num < 1)
			{
				throw invalid_argument("Invalid vector hash function value.");
			}
		}
	}

	// If values are not given
	if (this->cluster_num == -1)
	{
		while( this->cluster_num <= 0)
		{
			cout << "Insert cluster number:" << endl;
			cin >> this->cluster_num;
		}
	}
	if (this->vector_hfunc_num == -1)
	{
		while( this->vector_hfunc_num <= 0)
		{
			cout << "Insert vector's hash functions number:" << endl;
			cin >> this->vector_hfunc_num;
		}
	}

	#if DEBUG
		cout << "Classification out" << endl;
	#endif
}

Classification::~Classification()
{
	#if DEBUG
		cout << "~Classification in" << endl;
	#endif

	delete_vector<Clustering>(&clusterings);

	#if DEBUG
		cout << "~Classification out" << endl;
	#endif
}

// Class Classification_Points functions
Classification_Points::Classification_Points(string input_file, string config, short int flag): Classification(config)
{
	#if DEBUG
		cout << "Classification_Points in" << endl;
	#endif

	// Read points
	if ( !read(input_file, &(this->data)) )
		throw;

	// Create LSH ANN
	if (flag == 0 || ((flag >> 1) % 2 == 1))
	{
		try
		{
			this->lsh = new LSH(&(this->data), this->vector_htables_num, this->vector_hfunc_num, this->data.at(0)->get_dimension());
		}
		catch(bad_alloc&)
		{
			cerr << "main: No memory_available" << endl;
			throw;
		}
	}

	// Initialize clusterings
	if (flag == 0)
	{
		for (int i = 0; i < 8; i++)
		{
			try
			{
				this->clusterings.push_back(new Point_Clustering(i, this->cluster_num, &(this->data), this->lsh));

			}
			catch (std::bad_alloc & ba)
			{
				cerr << "Problem in Clustering initialization." << endl;
				throw;
			}
		}

	}
	else
	{
		try
		{
			this->clusterings.push_back(new Point_Clustering(flag, this->cluster_num, &(this->data), this->lsh));
		}
		catch (std::bad_alloc & ba)
		{
			cerr << "Problem in Clustering initialization." << endl;
			throw;
		}
	}

	#if DEBUG
		cout << "Classification_Points out" << endl;
	#endif

}

Classification_Points::~Classification_Points()
{
	#if DEBUG
		cout << "~Classification_Points in" << endl;
	#endif

	this->data.clear();	// No need to delete points. They are going to be deleted by lsh

	delete this->lsh;

	#if DEBUG
		cout << "~Classification_Points out" << endl;
	#endif
}

// Class Classification_Curves functions
Classification_Curves::Classification_Curves(string input_file, string config, short int flag): Classification(config)
{
	#if DEBUG
		cout << "Classification_Curves in" << endl;
	#endif

	this->grid_lsh = NULL;

	ifstream data;

	data.open(input_file);
	int i = 0;
	string line;
	int min_d = -1;	// Min curve dimension
	int max_d = -1;	// Max curve dimension
	if (data.is_open())
	{
		getline(data, line);	// Ignore type identification line

		// Get line. Each line is a curve
		while ( getline (data, line) && line.length() > 0 )
		{
			int pos1, pos2;
			string sub;

			// Find curve ID
			pos2 = line.find("\t");
			string id = line.substr(0, pos2);

			// Find curve size
			line = line.substr(pos2+1);
			pos2 = line.find("\t");
			sub = line.substr(0, pos2);

			// Add curve to data vector
			this->data.push_back(new Curve(id));

			if (line.empty())
			{
				delete_vector <Curve>(&(this->data));
				throw invalid_argument("Invalid form of curve");
			}

			// Find curve points
			double x,y;
			double last_x, last_y;
			bool first = true;
			while (!line.empty())
			{
				// Find coordinate x
				pos1 = line.find("(");
				if (pos1 < 0) break;
				pos2 = line.find(",");
				sub = line.substr(pos1 + 1, pos2 - pos1 - 1);
				x = stod(sub);

				// Move line
				line = line.substr(pos2 + 2);

				// Find coordinate y
				pos2 = line.find(")");
				sub = line.substr(0, pos2);
				y = stod(sub);

				// Move line
				line = line.substr(pos2 + 1);

				// Add point in curve if its different than the last one
				if (first)
				{
					this->data.at(i)->add_point(x, y);
					last_x = x;
					last_y = y;
					first = false;
				}
				else if (last_x != x || last_y != y)
				{
					this->data.at(i)->add_point(x, y);
					last_x = x;
					last_y = y;
				}
			}

			if (min_d == -1 && max_d == -1)
			{
				min_d = this->data.at(i)->get_length();
				max_d = this->data.at(i)->get_length();
			}
			else if (min_d > this->data.at(i)->get_length())
			{
				min_d = this->data.at(i)->get_length();
			}
			else if (max_d < this->data.at(i)->get_length())
			{
				max_d = this->data.at(i)->get_length();
			}

			i++;
		}

		data.close();
	}

	// Create grid_lsh structures
	if (flag == 0 || ((flag >> 1) % 2 == 1))
	{
		try
		{
			this->grid_lsh = new Grid_LSH(&(this->data), this->vector_htables_num, this->vector_hfunc_num, max_d, min_d);
		}
		catch(bad_alloc&)
		{
			cerr << "main: No memory available" << endl;
			throw;
		}
	}

	// Initialize clusterings
	if (flag == 0)
	{
		for (int i = 0; i < 8; i++)
		{
			try
			{
				this->clusterings.push_back(new Curve_Clustering(i, this->cluster_num, &(this->data), this->grid_lsh, min_d, max_d));
			}
			catch (std::bad_alloc & ba)
			{
				cerr << "Problem in Clustering initialization." << endl;
				throw;
			}
		}

	}
	else
	{
		try
		{
			this->clusterings.push_back(new Curve_Clustering(flag, this->cluster_num, &(this->data), this->grid_lsh, min_d, max_d));
		}
		catch (std::bad_alloc & ba)
		{
			cerr << "Problem in Clustering initialization." << endl;
			throw;
		}
	}

	#if DEBUG
		cout << "Classification_Curves out" << endl;
	#endif
}

Classification_Curves::~Classification_Curves()
{
	#if DEBUG
		cout << "~Classification_Curves in" << endl;
	#endif

	this->data.clear();	// No need to delete curves. They are going to be deleted by lsh

	delete this->grid_lsh;

	#if DEBUG
		cout << "~Classification_Curves in" << endl;
	#endif
}

// Class Clustering functions
Clustering::Clustering(int cluster_num)
{
	#if DEBUG
		cout << "Clustering in" << endl;
	#endif

	// The centers belong in the dataset
	for (int i = 0; i < cluster_num; ++i)
	{
		this->mean_centers.push_back(false);
	}

	#if DEBUG
		cout << "Clustering out" << endl;
	#endif
}

Clustering::~Clustering()
{
	#if DEBUG
		cout << "~Clustering in" << endl;
	#endif

	this->mean_centers.clear();
	this->clusters.clear();

	#if DEBUG
		cout << "~Clustering out" << endl;
	#endif
}

template<typename vector_type>
void Clustering::initialization1(unsigned int cluster_num, vector <vector_type *> *centers, vector<vector_type*>* data)
{
	#if DEBUG
		cout << "Clustering::initialization1 in" << endl;
	#endif

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, data->size()-1);

    int new_center;
    bool append;

	while(centers->size() < cluster_num)
	{
		new_center = dis(gen);
		append = true;

		// Check if its already a center or if there is a center with distance 0 from the new center
		for (unsigned int i = 0; i < centers->size(); i++)
		{
			if (this->distance(centers->at(i), data->at(new_center)) == 0)
			{
				append = false;
				break;
			}
		}

		// If the new center is approved, insert in centers' set
		if (append)
		{
			centers->push_back(data->at(new_center));
		}
	}

	#if DEBUG
		cout << "Clustering::initialization1 out" << endl;
	#endif
}

template <typename lsh_type, typename vector_type>
void Clustering::assignment2(lsh_type * lsh, vector <vector_type *> * data, vector <vector_type *> *centers)
{
	#if DEBUG
		cout << "Clustering::assignment2 in" << endl;
	#endif

	if (centers->size() == 0)
	{
		return;
	}

	// Use lsh to find nearest (to each center) Points/Curves
	unordered_map <vector_type *, int> clusters_temp;
	for (unsigned int cluster = 0; cluster < centers->size(); cluster++)
	{
	 	vector <vector_type *> * bucket = lsh->get_bucket(centers->at(cluster));

	 	if (bucket == NULL)
	 	{
	 		continue;
	 	}

	 	for (unsigned int j = 0; j < bucket->size(); j++)
	 	{
	 		if (cluster == 0)
	 		{
	 			clusters_temp.insert(make_pair(bucket->at(j), cluster));
	 		}
	 		else if (clusters_temp.count(bucket->at(j)))
	 		{
	 			if (this->distance(bucket->at(j), centers->at(clusters_temp[bucket->at(j)])) > this->distance(bucket->at(j), centers->at(cluster)))
	 			{
	 				clusters_temp[bucket->at(j)] = cluster;
	 			}
	 		}
	 		else
	 		{
	 			clusters_temp.insert(make_pair(bucket->at(j), cluster));
	 		}
	 	}

	 	delete bucket;
	}

	#if DEBUG
		cout << "Clustering::assignment2 Finished lsh" << endl;
	#endif

	this->clusters.clear();

	// Add the rest of Points/Curves to clusters
	for (unsigned int i = 0; i < data->size(); i++)
	{
		if (clusters_temp.count(data->at(i)))
		{
			this->clusters.insert(make_pair(clusters_temp[data->at(i)], i));
		}
		else
		{
			double min_dist = this->distance(data->at(i), centers->at(0));
			int min_cluster = 0;
			for (unsigned int cluster = 1; cluster < centers->size() ; cluster++)
			{
				double new_dist = this->distance(data->at(i), centers->at(cluster));
				if (new_dist < min_dist)
				{
					min_dist = new_dist;
					min_cluster = cluster;
				}
			}
			this->clusters.insert(make_pair(min_cluster, i));
		}
	}

	clusters_temp.clear();

	#if DEBUG
		cout << "Clustering::assignment2 out" << endl;
	#endif
}

// Class Point_Clustering functions
Point_Clustering::Point_Clustering(short int flag, int cluster_num, vector<Point*>* data, LSH * lsh) : Clustering(cluster_num)
{
	#if DEBUG
		cout << "Point_clustering in" << endl;
	#endif

	this->flag = flag;

	if (this->flag % 2 == 0)	// xx0 == Initialization 1
	{
		this->initialization1(cluster_num, &(this->centers), data);
	}
	else	// xx1 == Initialization 2
	{
		// this->initialization2(cluster_num, data);
	}

	int reps = 20;	// Maximum number of updates
	if (this->flag < 4)	// 0xx == Update 1
	{
		// TODO
	}
	if (this->flag >= 4)	// 1xx == Update 2
	{
		do
		{
			if ((this->flag >> 1) % 2 == 0)	// x0x == Assign 1
			{

			}
			else	// x1x == Assign 2
			{
				assignment2(lsh, data, &(this->centers));
			}

			#if DEBUG
			for (unsigned int i = 0; i < this->centers.size(); ++i)
			{
				cout << this->centers.at(i)->get_id() << endl;
				cout << this->clusters.count(i) << endl;
			}
			getchar();
			#endif
		}
		while (reps-- && update2(data));
	}

	#if DEBUG
		cout << "Point_Clustering out" << endl;
	#endif
}

Point_Clustering::~Point_Clustering()
{
	for (unsigned int i = 0; i < this->centers.size(); i++)
	{
		if (this->mean_centers.at(i))
		{
			delete this->centers.at(i);
		}
	}

	this->centers.clear();
}

void Point_Clustering::initialization2(int cluster_num, vector<Point*>* data){}

void Point_Clustering::assignment1(vector<Point*>* data){

	// Erase the previous clusters
	// this->clusters.clear();

	// for (unsigned int i = 0; i < (*data).size(); ++i)
	// {
	// 	double min_dist = -1, cur_dist;
	// 	int min_c;
	// 	if ( this->centers.find(i) != this->centers.end() )		// If the point is not a center
	// 	{
	// 		set <int, greater <int> > :: iterator itr;
	// 		for (itr = this->centers.begin(); itr != this->centers.end(); ++itr)	// For every center
	// 		{
	// 			if ( min_dist == -1 )	// If this is the first distance calculated
	// 			{
	// 				min_dist = manhattan_dist((*data)[*itr], (*data)[i]);
	// 				min_c = *itr;
	// 			}
	// 			else
	// 			{
	// 				cur_dist = manhattan_dist((*data)[*itr], (*data)[i]);	// Calculate distance between
	// 				if ( min_dist < cur_dist )	// If the current center is closer
	// 				{
	// 					min_dist = cur_dist;
	// 					min_c = *itr;
	// 				}
	// 			}
	// 		}
	// 	}
	// 	this->clusters.insert(pair<int, int> (min_c, i));
	// }
}

bool Point_Clustering::update2(vector<Point*>* data)
{
	#if DEBUG
		cout << "Point_Clustering::update2 in" << endl;
	#endif

	vector <Point *> new_centers;
	vector <bool> new_mean_centers;

	if (this->centers.size() == 0)
	{
		return false;
	}

	// For each cluster find mean vector
	for (unsigned int i = 0; i < this->centers.size(); i++)
	{
		int cluster_size = this->clusters.count(i);	// If cluster is empty assign a new center found randomly

		// Assign random center
		if (cluster_size == 0)
		{
			/* TODO We are using previous center. Not a random one */
			new_centers.push_back(NULL);
			new_mean_centers.push_back(false);
			continue;
		}

		// Create mean vector
		vector <double> mean_vector;
		for (int j = 0; j < data->at(i)->get_dimension(); j++)
		{
			mean_vector.push_back(0.0);
		}

		pair <multimap<int, int>::iterator, multimap<int, int>::iterator> range = this->clusters.equal_range(i);
		for (multimap<int, int>::iterator it = range.first; it != range.second; ++it)
		{
			for (unsigned int j = 0; j < mean_vector.size(); j++)
			{
				mean_vector.at(j) += (*(data->at(it->second)))[j];
			}
		}

		// Divide with cluster size
		for (unsigned int j = 0; j < mean_vector.size(); j++)
		{
			mean_vector.at(j) = mean_vector.at(j)/cluster_size;
		}
		new_centers.push_back(new Point("none", &mean_vector));
		new_mean_centers.push_back(true);

	}

	// TODO  synartisi annas
	// Compare old centers with new centers
	// If you have to continue return true else false

	// bool changed = synarthsh annas

	bool changed = false;
	for (unsigned int i = 0; i < this->centers.size(); i++)
	{
		if (new_centers.at(i) != NULL && this->distance(this->centers.at(i), new_centers.at(i)))
		{
			changed = true;
		}
	}

	// Replace the centers if needed
	if (changed)
	{
		for (unsigned int i = 0; i < this->centers.size(); ++i)
		{
			// Delete non dataset centers if they have been changed
			if (this->clusters.count(i) && this->mean_centers.at(i))
			{
				delete this->centers.at(i);
				this->centers.at(i) = new_centers.at(i);
			}
			else if (new_centers.at(i) != NULL)
			{
				this->centers.at(i) = new_centers.at(i);
				this->mean_centers.at(i) = new_mean_centers.at(i);
			}
		}
	}
	else
	{
		for (unsigned int i = 0; i < new_centers.size(); ++i)
		{
			delete new_centers.at(i);
		}
		new_centers.clear();
		new_mean_centers.clear();
	}

	#if DEBUG
		cout << "Point_Clustering::update2 out" << endl;
	#endif

	return changed;
}

double Point_Clustering::binary_search(vector<double>* P, double x)
{
	int l = 0, h = (*P).size()-1, m;
	double middle_value;

	// Corner case
	if ( x == (*P)[h] )
		return h;

	// Binary search
	while ( l < h )
	{

		m = (h+l)/2;
		middle_value = (*P)[m];
		if ( middle_value == x )
			return m;
		else if ( x < middle_value )
		{
			if ( m > 0 && x > (*P)[m-1]) 	// If this is not the first element and P[m-1] < x < P[m]
                return m;
			h = m;		// Choose the left half
		}
		else
		{
			if ( m < (int)((*P).size()-1) && x < (*P)[m+1] )	// If this is not the last element and P[m] < x < P[m+1]
                return m+1;
			l = m + 1;	// Choose the right half
		}
	}

	return m;
}

double Point_Clustering::distance(Point *c1, Point *c2)
{
	return manhattan_dist(c1,c2);
}

// Class Curve_Clustering functions
Curve_Clustering::Curve_Clustering(short int flag, int cluster_num, vector<Curve*>* data, Grid_LSH * grid_lsh, int min_d, int max_d) : Clustering(cluster_num)
{
	#if DEBUG
		cout << "curves_clustering in" << endl;
	#endif

	this->flag = flag;

	if (this->flag % 2 == 0)	// xx0 == Initialization 1
	{
		this->initialization1(cluster_num, &(this->centers), data);
	}
	else	// xx1 == Initialization 2
	{
		// this->initialization2(cluster_num, data);
	}

	int reps = 20;	// Maximum number of updates
	if (this->flag < 4)	// 0xx == Update 1
	{
		// TODO
	}
	if (this->flag >= 4)	// 1xx == Update 2
	{
		do
		{
			if ((this->flag >> 1) % 2 == 0)	// x0x == Assign 1
			{

			}
			else	// x1x == Assign 2
			{
				assignment2(grid_lsh, data, &(this->centers));
			}
		}
		while (reps-- && update2(data));
	}

	#if DEBUG
		cout << "Curve_Clustering out" << endl;
	#endif
}

Curve_Clustering::~Curve_Clustering()
{
	#if DEBUG
		cout << "Curve_Clustering in" << endl;
	#endif

	for (unsigned int i = 0; i < this->centers.size(); i++)
	{
		if (this->mean_centers.at(i))
		{
			delete this->centers.at(i);
		}
	}

	this->centers.clear();

	#if DEBUG
		cout << "Curve_Clustering out" << endl;
	#endif
}

bool Curve_Clustering::update2(vector<Curve*>* data)
{
	#if DEBUG
		cout << "Curve_Clustering::update2 in" << endl;
	#endif

	vector <Curve *> new_centers;

	if (this->centers.size() == 0)
	{
		return false;
	}

	for (unsigned int cluster = 0; cluster < this->centers.size(); cluster++)
	{
		int cluster_size = this->clusters.count(cluster);	// If cluster is empty assign a new center found randomly

		// Assign random center
		if (cluster_size == 0)
		{
			/* TODO */
			new_centers.push_back(NULL);
			continue;
		}


		pair <multimap <int ,int>::iterator, multimap <int, int>::iterator> range;
    	range = this->clusters.equal_range(cluster);

		// Calculate mean size of curve in cluster (λ)
		int mean_dim = 0;
		for (multimap<int,int>::iterator it=range.first; it!=range.second; ++it)
		{
			mean_dim += data->at(it->second)->get_length();
		}

		mean_dim = (int) floor((double) mean_dim/cluster_size);

		int counter = 0;	// How many curves have size >= mean_dim
		for (multimap<int,int>::iterator it=range.first; it!=range.second; ++it)
		{
			if (data->at(it->second)->get_length() >= mean_dim)
			{
				counter++;
			}
		}

		// Find random curve of size >= mean_dim
		random_device rd;  //Will be used to obtain a seed for the random number engine
    	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    	uniform_int_distribution<int> dis(0, counter-1);
		int random_curve = dis(gen);
		int random_curve_pos;
		for (multimap<int,int>::iterator it=range.first; it!=range.second; ++it)
		{
			if (data->at(it->second)->get_length() >= mean_dim)
			{
				random_curve--;
				random_curve_pos = it->second;
			}

			if (random_curve == 0) break;
		}

		// Create the curve for the algorithm
		Curve * mean_curve = NULL;
		Curve * new_mean_curve = new Curve("none");
		int start_point = 0;
		// If random curve is bigger than mean_dim, choose a random section of size mean_dim
		if (data->at(random_curve_pos)->get_length() > mean_dim)
		{
			start_point = rand() % (data->at(random_curve_pos)->get_length() - mean_dim);
		}

		for (int i = 0; i < mean_dim; ++i)
		{
			new_mean_curve->add_point((*(data->at(random_curve_pos)))[start_point+i].first, (*(data->at(random_curve_pos)))[start_point+i].second);
		}

		do
		{
			delete mean_curve;
			mean_curve = new_mean_curve;
			vector <pair <double, double>> * A = new vector <pair <double, double>>	 [mean_dim];
			for (multimap<int,int>::iterator it=range.first; it!=range.second; ++it)
			{
				vector <pair <int, int>> opt_trav;
				DTW_distance(mean_curve, data->at(it->second), &opt_trav);
				for (unsigned int i = 0; i < opt_trav.size(); ++i)
				{
					A[opt_trav.at(i).first].push_back((*(data->at(it->second)))[opt_trav.at(i).second]);
				}
			}

			// Calculate mean and create new Curve (C)
			new_mean_curve = new Curve("");
			for (int i = 0; i < mean_dim; i++)
			{
				double x = 0, y = 0;
				for (unsigned int j = 0; j < A[i].size(); j++)
				{
				  	x += A[i].at(j).first;
				  	y += A[i].at(j).second;
				}
				if (A[i].size() > 0)
				{
					x = x/A[i].size();
					y = x/A[i].size();
				}
				new_mean_curve->add_point(x, y);
			}

		} while (DTW_distance(mean_curve, new_mean_curve) == 0);

		delete new_mean_curve;

		// Add new center to dataset
		new_centers.push_back(mean_curve);
	}

	// TODO  synartisi annas
	// Compare old centers with new centers
	// If you have to continue return true else false

	// bool changed = synarthsh annas

	bool changed = false;
	for (unsigned int i = 0; i < this->centers.size(); i++)
	{
		if (new_centers.at(i) != NULL && this->distance(this->centers.at(i), new_centers.at(i)))
		{
			changed = true;
		}
	}

	// Replace the centers if needed
	if (changed)
	{
		for (unsigned int i = 0; i < this->centers.size(); ++i)
		{
			// Delete non dataset centers if they have been changed
			if (this->clusters.count(i) && this->mean_centers.at(i))
			{
				delete this->centers.at(i);
				this->centers.at(i) = new_centers.at(i);
			}
			else if (this->clusters.count(i))
			{
				this->centers.at(i) = new_centers.at(i);
				this->mean_centers.at(i) = true;
			}
		}
	}
	else
	{
		for (unsigned int i = 0; i < new_centers.size(); ++i)
		{
			delete new_centers.at(i);
		}
		new_centers.clear();
	}

	#if DEBUG
		cout << "Curve_Clustering::update2 out" << endl;
	#endif

	// return changed;
	return true;
}

double Curve_Clustering::distance(Curve *c1, Curve *c2)
{
	return DTW_distance(c1,c2);
}
