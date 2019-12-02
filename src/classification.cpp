#include "../include/classification.hpp"

#include <fstream>
#include <random>

// Class Classification functions
Classification::Classification(string config)
{
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

}

Classification::~Classification()
{
	delete_vector<Clustering>(&clusterings);
}

// void Classification::train()
// {
// }

// Class Classification_Points functions
Classification_Points::Classification_Points(string input_file, string config, short int flag): Classification(config)
{
	// Read points
	if ( !read(input_file, &(this->data)) )
		throw;

	// Initialize clusterings
	if (flag == 0)
	{	
		for (int i = 0; i < 8; i++)
		{
			try
			{
				this->clusterings.push_back(new Point_Clustering(i, this->cluster_num, &(this->data)));

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
			this->clusterings.push_back(new Point_Clustering(flag, this->cluster_num, &(this->data)));
		}
		catch (std::bad_alloc & ba)
		{
			cerr << "Problem in Clustering initialization." << endl;
			throw;
		}
	}

}

Classification_Points::~Classification_Points()
{
	delete_vector<Point>(&(this->data));
}

// Class Classification_Curves functions
Classification_Curves::Classification_Curves(string input_file, string config, short int flag): Classification(config)
{
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
	try
	{
		this->grid_lsh = new Grid_LSH(&(this->data), this->vector_htables_num, this->vector_hfunc_num, max_d, min_d);
	}
	catch(bad_alloc&)
	{
		cerr << "main: No memory available" << endl;
		throw;
	}

	// Initialize clusterings
	if (flag == 0)
	{
		for (int i = 0; i < 8; i++)
		{
			try
			{
				this->clusterings.push_back(new Curve_Clustering(i, this->cluster_num, &(this->data), min_d, max_d));
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
			this->clusterings.push_back(new Curve_Clustering(flag, this->cluster_num, &(this->data), min_d, max_d));
		}
		catch (std::bad_alloc & ba)
		{
			cerr << "Problem in Clustering initialization." << endl;
			throw;
		}
	}
}

Classification_Curves::~Classification_Curves()
{
	delete_vector<Curve>(&(this->data));
}

// Class Clustering functions
Clustering::~Clustering()
{
	this->clusters.clear();
}

template<typename vector_type, typename Function>
void Clustering::initialization1(unsigned int cluster_num, vector <vector_type *> *centers, vector<vector_type*>* data, Function dist_function)
{
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
	// this->distance(centers.at(i).second), data->at(new_center));
			if ((dist_function(centers->at(i), data->at(new_center)) == 0))
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
}

// Class Point_Clustering functions
Point_Clustering::Point_Clustering(short int flag, int cluster_num, vector<Point*>* data) /*:Clustering(flag, cluster_num, data)*/{

	this->flag = flag;

	if (this->flag % 2 == 0)	// xx0 == Initialization 1
	{
		this->initialization1(cluster_num, &(this->centers), data, manhattan_dist);
	}
	else	// xx1 == Initialization 2
	{
		// this->initialization2(cluster_num, data);
	}

	if (this->flag < 4)	// 0xx == Update 1
	{
		// TODO
	}
	if (this->flag >= 4)	// 1xx == Update 2
	{
		while (update2(data))
		{
			// x0x == Assign 1
			// x1x == Assign 2
		}
	}
}

Point_Clustering::~Point_Clustering()
{
	while (!this->centers.empty())
	{
		if (this->centers.back()->get_id().compare("none") == 0)
		{
			delete this->centers.back();
		}
		this->centers.pop_back();
	}
}

// Class Curve_Clustering functions
Curve_Clustering::Curve_Clustering(short int flag, int cluster_num, vector<Curve*>* data, int min_d, int max_d) /*:Clustering(flag, cluster_num, data)*/
{
cout << "curves_clustering in" << endl;
	this->flag = flag;

	if (this->flag % 2 == 0)	// xx0 == Initialization 1
	{
		this->initialization1(cluster_num, &(this->centers), data, DTW);
	}
	else	// xx1 == Initialization 2
	{
		// this->initialization2(cluster_num, data);
	}

	
	if (this->flag < 4)	// 0xx == Update 1
	{
		// TODO
	}
	if (this->flag >= 4)	// 1xx == Update 2
	{
		while (update2(data, min_d, max_d))
		{
			// x0x == Assign 1
			// x1x == Assign 2
		}
	}
	
cout << "Curve_Clustering out" << endl;
}

Curve_Clustering::~Curve_Clustering()
{
	while (!this->centers.empty())
	{
		if (this->centers.back()->get_id().compare("none") == 0)
		{
			delete this->centers.back();
		}
		this->centers.pop_back();
	}
}

void Point_Clustering::initialization2(int cluster_num, vector<Point*>* data){

	// Δεν νομιζω πως χρειαζεται γιατι θα ειναι αδεια
// 	this->centers.clear();

// 	int data_size = data->size();

// // cout << "\ncluster_num is " << cluster_num << endl;

// 	random_device rd;  //Will be used to obtain a seed for the random number engine
//     mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//     uniform_int_distribution<> dis(0, data_size-1);

// 	this->centers.insert(dis(gen));

// 	for (int t = 1; t < cluster_num; ++t)	// Until all centers have been chosen
// 	{
// 		double sum = 0;
// 		vector<double> P;

// 		double max_D_i = 0;
// 		for (int i = 0; i < data_size; ++i)	// For every point
// 		{
// 			if ( this->centers.find(i) == this->centers.end() )	// if the point is not a center
// 			{
// 				double D_i = this->min_dist(data, i);
// 				sum += D_i*D_i;
// 				P.push_back(sum);

// 				if ( max_D_i < D_i )
// 					max_D_i = D_i;
// 			}
// 		}

// 		// Normalize P
// 		for (unsigned int i = 0; i < P.size(); ++i)
// 		{
// 			P[i] = P[i]/max_D_i; 
// 		}

//     	uniform_real_distribution<> dis(0.0, P[P.size()-1]/1.0);

//     	double x = dis(gen);
// // cout << "x is " << x << " and p is " << P[P.size()-1] <<endl;
//    		int r = this->binary_search(&P, x);	// P(r-1) < x <= P(r)
// // if ( r!=0 )cout << "p(r-1) is " << P[r-1] << endl;
// // cout << "p(r) is " << P[r] << endl;
//    		int c_num = 0;

//    		// Find num of centers before the r_pos
// 		set <int, greater <int> > :: iterator itr = centers.begin();
//    		for (itr = this->centers.begin(); itr != this->centers.end(); ++itr)
//    		{
//    			if ( (*itr) <= r )
//    				c_num++;
//    		}

//    		while ( this->centers.find(r+c_num) != this->centers.end() ) c_num++;	// While the current point is a center
//     	this->centers.insert(r+c_num);
// 	}

// // cout << "centers are " << endl;
// // 	set <int, greater <int> > :: iterator itr = centers.begin();
// //    	for (itr = this->centers.begin(); itr != this->centers.end(); ++itr)
// //    	{
// //    		cout << *itr << endl;
// //    	}
}

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
	vector <pair <int, Point *>> new_center;

	// For each cluster find mean vector
	for (unsigned int i = 0; i < this->centers.size(); i++)
	{
		int cluster_size = 0;	// If cluster is empty assign a new center found randomly

		// Create mean vector
		vector <double> mean_vector;
		for (int i = 0; i < data->at(i)->get_dimension(); ++i)
		{
			mean_vector.push_back(0.0);
		}

		pair <multimap<int, int>::iterator, multimap<int, int>::iterator> range = this->clusters.equal_range(i);
		for (multimap<int, int>::iterator it = range.first; it != range.second; ++it)
		{
			cluster_size++;
			for (unsigned int i = 0; i < mean_vector.size(); ++i)
			{
				mean_vector.at(i) += (*data->at(it->second))[i];			
			}
		}

		// Assign random center
		if (cluster_size == 0)
		{
			/* TODO */
		}
		else
		{
			// Divide with cluster size
			for (unsigned int i = 0; i < mean_vector.size(); ++i)
			{
				mean_vector.at(i) = mean_vector.at(i)/cluster_size;
			}

			new_center.push_back(make_pair(i, new Point("none", &mean_vector)));
		}
	}

	// TODO
	// Compare old centers with new centers
	// If you have to continue return true else false

	return true;
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

// double Point_Clustering::min_dist(vector<Point*>* data, int pos)
// {
// 	Point* point = (*data)[pos];
// 	set <int, greater <int> > :: iterator itr;

// 	double min, cur_dist;
// 	for (itr = this->centers.begin(); itr != this->centers.end(); ++itr)
// 	{
// 		if ( itr == this->centers.begin() )	// If this is the first center
// 		{
// 			min = manhattan_dist(point, (*data)[*itr]);
// 		}
// 		else
// 		{
// 			cur_dist = manhattan_dist(point, (*data)[*itr]);
// 			if ( cur_dist < min )	// If the distance from the current center is smaller
// 				min = cur_dist;	
// 		}

// 	}

// 	return min;
// }


double Point_Clustering::distance(Point *c1, Point *c2)
{
	return manhattan_dist(c1,c2);
}

double Curve_Clustering::distance(Curve *c1, Curve *c2)
{
	return DTW_distance(c1,c2);
}