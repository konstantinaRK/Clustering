#include <iostream>
#include <string>
#include "../include/classification.hpp"
#include "../include/utilities.hpp"

using namespace std;

int main(int argc, char * argv [])
{
	string output_file;
	bool complete;

	Classification * classification = dataHandling(argc, argv, &output_file, &complete);
	if (classification == NULL)
	{
		return 1;
	}
cout << "lalaal" << endl;
	delete classification;
cout << "yoyoyo" << endl;

	return 0;
}
