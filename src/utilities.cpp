#include "../include/utilities.hpp"

#include <iostream>
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
