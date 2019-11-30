#include "../include/utilities.hpp"

#include <string.h>
#include <fstream>

using namespace std;

// Point
Point::Point(string id, vector<double>* X){

	this->id = id;
	for (unsigned int i=0; i<(*X).size(); i++){
		this->X.push_back((*X)[i]);
	}
}

Point::~Point(){
	this->X.clear();
}

int Point::get_dimension(void){
	return this->X.size();
}

double Point::operator[](unsigned int i){
	 if ( i>=0 && i<this->X.size() )
	 	return this->X[i];
	 else
	 	return -1; 
}

// Curve
Curve::Curve(string id){
	this->curve_id = id;
}

Curve::~Curve()
{
	this->points.clear();
}

void Curve::add_point(double x, double y){

	this->points.push_back(make_pair(x, y));
}

int Curve::get_length(void){

	return this->points.size();
}

pair<double, double> Curve::operator[](int i){

	return this->points[i];
}


void Curve::clear(void){

	this->points.clear();
}

// Functions
Classification * dataHandling(int argc, char * argv[], string * output_file, bool * complete)
{
	string data_file = "";
	string config = "";
	short int flag = 0;

	*output_file = "";
	*complete = false;

	// Check if the arguments are right and initialize values of variables (input_file, argument_file, output_file, complete print)
	if (argc < 7)
	{
		cerr << "Program must be executed as:\n./cluster –i <input file> –c <configuration file> -o <output file> -complete (optional)" << endl;
		return NULL;
	}

	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "-i") == 0 && i+1)
		{
			if (i+1 > argc)
			{
				cerr << "Program must be executed as:\n./cluster –i <input file> –c <configuration file> -o <output file> -complete (optional)" << endl;
				return NULL;
			}
			data_file = argv[i+1];
			i++;
		}
		else if (strcmp(argv[i], "-c") == 0)
		{
			if (i+1 > argc)
			{
				cerr << "Program must be executed as:\n./cluster –i <input file> –c <configuration file> -o <output file> -complete (optional)" << endl;
				return NULL;
			}
			config = argv[i+1];
			i++;
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			if (i+1 > argc)
			{
				cerr << "Program must be executed as:\n./cluster –i <input file> –c <configuration file> -o <output file> -complete (optional)" << endl;
				return NULL;
			}
			*output_file = argv[i+1];
			i++;
		}
		else if (strcmp(argv[i], "-complete") == 0)
		{
			*complete = true;
		}
		else if(strcmp(argv[i], "-interactive") == 0)
		{
			if (i+1 > argc)
			{
				cerr << "Program must be executed as:\n./cluster –i <input file> –c <configuration file> -o <output file> -complete (optional)" << endl;
				return NULL;
			}
			flag = stoi(argv[i+1]);
			if (flag < 0 || !(flag == 0 || flag == 111 || flag == 121 || flag == 112 || flag == 122 || flag == 211 || flag == 212 || flag == 221 || flag == 222))
			{
				cerr << "Flag must describe the combination of initialization|assignement|update" << endl;
				return NULL;
			}
			i++;
		}
	}

	if (data_file.empty() || config.empty() || (*output_file).empty())
	{
		cerr << "Program must be executed as:\n./cluster –i <input file> –c <configuration file> -o <output file> -complete (optional)" << endl;
		return NULL;
	}

	// Read input data from file
	ifstream data;
	data.open(data_file);
	if (!data.is_open())
	{
		cerr << "Invalid input file" << endl;
		return NULL;
	}

	Classification * return_value = NULL;
	string line;
	getline(data, line);
	data.close();
	if (!line.compare(0, 7, "vectors"))
	{
		try 
		{
			return_value = new Classification_Points(data_file, config, flag);
		}
		catch (const std::invalid_argument& ia){}
	}
	else if(!line.compare(0, 6, "curves"))
	{
		try 
		{
			return_value = new Classification_Curves(data_file, config, flag);
		}
		catch (const std::invalid_argument& ia){}
	}

	return return_value;
}

double DTW_distance(Curve* x1, Curve* x2, vector<pair<int, int>>* opt_trav){

	int m1 = x1->get_length();
	int m2 = x2->get_length();

	pair<double, string>** C;
	try{
		C = new pair<double, string>*[m1];
		for (int i = 0; i < m1; ++i)
		{
			C[i] = new pair<double, string>[m2];
		}
	}
	catch(bad_alloc&){
		cerr << "DTW_distance: No memory available" << endl;
		return -1;
	}

	C[0][0].first = eucl_dist((*x1)[0], (*x2)[0]);
	C[0][0].second = "";

	// Fill the first row
	for (int j = 1; j < m2; ++j)
	{
		C[0][j].first = C[0][j-1].first + eucl_dist((*x1)[0], (*x2)[j]);
		C[0][j].second = "left";
	}

	// Fill the first column
	for (int i = 1; i < m1; ++i){
		C[i][0].first = C[i-1][0].first + eucl_dist((*x1)[i], (*x2)[0]);
		C[i][0].second = "up";
	}

	// Fill the rest of the matrix
	string direction = "";
	for (int i = 1; i < m1; ++i)
	{
		for (int j = 1; j < m2; ++j)
		{
			C[i][j].first = min(C[i-1][j].first, C[i-1][j-1].first, C[i][j-1].first, &direction) + eucl_dist((*x1)[i], (*x2)[j]);
			C[i][j].second = direction;
		}
	}

	double distance = C[m1-1][m2-1].first;

	if ( opt_trav!= NULL )
	{
		// Find optimal traversal
		string direction;
		pair<int, int> pos;
		pos.first = m1;
		pos.second = m2;
		while ( pos.first!=0 && pos.second!=0 )
		{
			(*opt_trav).push_back(pos);
			direction = C[pos.first][pos.second].second;
			if ( direction.compare("left") == 0 )
				pos.first--;
			else if ( direction.compare("diag") == 0 )
			{
				pos.first--;
				pos.second--;
			}
			else // direction = up
				pos.second--;
		}
		reverse((*opt_trav).begin(),(*opt_trav).end());
	}

	// Free table
	for (int i = 0; i < m1; ++i)
	{
		delete [] C[i];
	}
	delete [] C;

	return distance;

}


double min(double x, double y, double z, string* direction){

	double min = x;
	if ( y < min && z >= y )
	{
		(*direction) = "diag";
		return y;
	}
	else if ( z < min )
	{
		(*direction) = "up";
		return z;
	}
	(*direction) = "left";
	return min;
}


double eucl_dist(pair<double, double> x, pair<double, double> y){

	double x1 = x.first - y.first;
	double x2 = x.second - y.second;
	return (sqrt((x1*x1)+(x2*x2)));
}

double manhattan_dist(Point* x, Point* y){

	double distance = 0;
	int xi, yi;
	double d;
	for (int i = 0; i < x->get_dimension(); ++i)
	{
		xi = (*x)[i];
		yi = (*y)[i];

		d = xi - yi;
		if ( d >= 0 )
			distance += d;
		else
			distance += yi - xi; 
	}

	return distance;
}

// πρεπει να τους αλλαξω μερος

// Read points from a file
bool read(string file_name, vector<Point*>* points){

	string line;
	int d;

  	ifstream myfile;
  	myfile.open(file_name);

  	if (myfile.is_open())
  	{
  		// Ignore the first line
  		getline(myfile, line);

  		getline(myfile, line);

  		//Read the first point and store d
    	if ( !point_proccessing(points, line) )
    		return false;
    	d = ((*points)[0])->get_dimension();

  		// Read rest points
    	while ( getline(myfile, line) ){
    	  if ( !point_proccessing(points, line, d) ){
    	  	return false;
    	  }
    	}
    	myfile.close();
  	}
  	else
  	{
  		cerr << "Unable to open file" << endl;
  		return false;
  	}
  	return true;
}

// Transform data to Point structure
bool point_proccessing(vector<Point*>* points, string p, int d){

	stringstream L;
	L << p;
	string coordinate;

	if ( !getline(L, coordinate, '\t') )
		return false;

	// Get the id
	string id = coordinate;

	// Get the coordinates
	vector<double> X;
	while( getline(L, coordinate, '\t') ){

		if ( coordinate.compare("\r") == 0 ) break;
		X.push_back( stod(coordinate) );
	}

	// If this is not the first point and the dimension is different in this point
	if ( d != -1 && (int)X.size() != d ){
		cerr << "Wrong input file" << endl;
		return false;
	}
	// Create a point
	Point* point;
	try{
		point = new Point(id, &X);
	}
	catch(std::bad_alloc&) {
		cerr << "No memory available" << endl;
		return false;   
	}
	X.clear();

	// Insert point in the list
	(*points).push_back(point);

	return true;
}


double DTW(Curve *x1, Curve * x2){

	return DTW_distance(x1, x2);
}