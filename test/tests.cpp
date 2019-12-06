// tests.cpp

#include "../include/classification.hpp"
#include "../include/utilities.hpp"
#include <gtest/gtest.h>

using namespace std;

TEST(manhattan_dist, NullPointers){

	ASSERT_EQ(-1, manhattan_dist(NULL, NULL));
	vector<double> vec{10, 20, 30, 40};
	Point* point = new Point("1", &vec);
	ASSERT_EQ(-1, manhattan_dist(NULL, point));
	ASSERT_EQ(-1, manhattan_dist(point, NULL));
	delete point;
}

TEST(manhattan_dist, SamePoints){

	vector<double> vec{1.5, 20, 3.7, 40};
	Point * point = new Point("1", &vec);
	ASSERT_EQ(0, manhattan_dist(point, point));
	delete point;
}

TEST(manhattan_dist, DifferentDimensions){

	vector<double> vec1{2, 5.8, 3, 40};
	vector<double> vec2{25, 2, 15, 4.5, 3.8};
	Point * point1 = new Point("1", &vec1), * point2 = new Point("2", &vec2);
	ASSERT_EQ(-1, manhattan_dist(point1, point2));
	delete point1;
	delete point2;
}

TEST(manhattan_dist, DifferentPoints){

	vector<double> vec1{2, 3};
	vector<double> vec2{3, 5};
	Point * point1 = new Point("1", &vec1), * point2 = new Point("2", &vec2);
	ASSERT_EQ(3, manhattan_dist(point1, point2));
	delete point1;
	delete point2;

	vector<double> vec3{2, 4.2, 8};
	vector<double> vec4{1, 5, 6.3};
	Point * point3 = new Point("1", &vec3), * point4 = new Point("2", &vec4);
	ASSERT_EQ(3.5, manhattan_dist(point3, point4));
	delete point1;
	delete point2;
}

TEST(DTW_distance, NullPointers){

	ASSERT_EQ(-1, DTW_distance(NULL, NULL));
	vector<pair<double, double>> vec{make_pair(1, 2), make_pair(2, 5), make_pair(3, 8) };
	Curve* curve = new Curve("1");
	for (unsigned int i = 0; i < vec.size(); ++i)
	{
		curve->add_point(vec[i].first, vec[i].second);
	}
	ASSERT_EQ(-1, DTW_distance(NULL, curve));
	ASSERT_EQ(-1, DTW_distance(curve, NULL));
	delete curve;
}

TEST(DTW_distance, SamePoints){

	vector<pair<double, double>> vec{make_pair(10, 20), make_pair(30, 40) };
	Curve* curve = new Curve("1");
	for (unsigned int i = 0; i < vec.size(); ++i)
	{
		curve->add_point(vec[i].first, vec[i].second);
	}
	ASSERT_EQ(0, DTW_distance(curve, curve));
	delete curve;
}

TEST(DTW_distance, DifferentPoints){

	vector<pair<double, double>> vec1{ make_pair(2, 3), make_pair(5, 8), make_pair(14, 0) };
	vector<pair<double, double>> vec2{ make_pair(3, 5), make_pair(6, 8) };
	Curve * curve1 = new Curve("1"), * curve2 = new Curve("2");
	for (unsigned int i = 0; i < vec1.size(); ++i)
	{
		curve1->add_point(vec1[i].first, vec1[i].second);
	}
	for (unsigned int i = 0; i < vec1.size(); ++i)
	{
		curve2->add_point(vec2[i].first, vec2[i].second);
	}
	ASSERT_EQ(17.2361, DTW_distance(curve1, curve2));
	delete curve1;
	delete curve2;

	vector<pair<double, double>> vec3{ make_pair(2, 3), make_pair(5.5, 10) };
	vector<pair<double, double>> vec4{ make_pair(1, 0.5), make_pair(2, 5) };
	curve1 = new Curve("1");
	curve2 = new Curve("2");
	for (unsigned int i = 0; i < vec3.size(); ++i)
	{
		curve1->add_point(vec3[i].first, vec3[i].second);
	}
	for (unsigned int i = 0; i < vec4.size(); ++i)
	{
		curve2->add_point(vec4[i].first, vec4[i].second);
	}
	ASSERT_EQ(8.79586, DTW_distance(curve1, curve2));
	delete curve1;
	delete curve2;
}

TEST(min, SameValues){

	string direction;
	ASSERT_EQ(1.2, min(1.2, 1.2, 1.2, &direction));
}

TEST(min, DifferentValues){

	string direction;
	ASSERT_EQ(1, min(1, 3, 2, &direction));
	ASSERT_EQ(-1, min(5, -1, 8, &direction));
	ASSERT_EQ(2, min(2, 2, 10, &direction));
}

// TEST(binary_search, NullPointer){

// 	ASSERT_EQ(-1, Clustering::binary_search(NULL, 2));
// }

// TEST(binary_search, CornerCases){

// 	vector<double> vec{0, 2, 5, 8, 15, 20};
// 	ASSERT_EQ(0, Clustering::binary_search(&vec, 0));
// 	ASSERT_EQ(5, Clustering::binary_search(&vec, 20));
// }

// TEST(binary_search, NormalCases){

// 	vector<double> vec{0, 2, 5, 8, 15, 20};
// 	ASSERT_EQ(2, Clustering::binary_search(&vec, 4));
// 	ASSERT_EQ(4, Clustering::binary_search(&vec, 14));
// }

int main(int argc, char **argv){

	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}