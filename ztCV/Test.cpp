#include "Types.h"
#include "Types.hpp"
#include "Mat.h"
#include "Mat.hpp"
#include <cstdio>
#include <cassert>
#include <iostream>
using namespace ztCV; 

static int main_ret = 0;
static int test_count = 0;
static int test_pass = 0;

#define EXPECT_EQ_BASE(equality, expect, actual, format)\
	do{\
		test_count++;\
		if(equality)\
			test_pass++;\
		else{\
			fprintf(stderr, "%s:%d: expect: " format ", actual: " format "\n", __FILE__, __LINE__, expect, actual);\
			main_ret = 1;\
		}\
	} while(0)

#define EXPECT_EQ_INT(expect, actual) EXPECT_EQ_BASE((expect) == (actual), expect, actual, "%d")
#define EXPECT_TRUE(actual) EXPECT_EQ_BASE((actual)!=0, "true", "false", "%s")
#define EXPECT_FALSE(actual) EXPECT_EQ_BASE((actual)==0, "false", "true", "%s")

#define EXPECT_EQ_POINT(expect, actual)\
	do{\
		assert(expect.x_ == actual.x_);\
		assert(expect.y_ == actual.y_);\
	} while(0)
	

static void test_point() {
	Point pt1(1,2);
	Point pt2(3, 4);
	Point pt3(pt1);
	EXPECT_EQ_INT(1, pt1.x_);
	EXPECT_EQ_INT(2, pt1.y_);
	EXPECT_EQ_INT(11, pt1.inner_product(pt2));
	EXPECT_EQ_INT(10, pt1.outer_product(pt2));
	EXPECT_EQ_POINT(pt1, pt3);
}

static void test_size() {
	Size sz(4, 8);
	EXPECT_EQ_INT(4, sz.width_);
	EXPECT_EQ_INT(8, sz.height_);
	EXPECT_EQ_INT(4 * 8, sz.area());
	EXPECT_FALSE(sz.empty());
}




#define TEST_MAT_ZEROS_BASE(rows, cols, type)\
	do{\
		Mat mat = Mat::zeros(rows, cols, type);\
		std::cout << "mat=" << std::endl;\
		std::cout << mat << std::endl;\
	} while(0)

#define TEST_MAT_ZEROS1(size, type) TEST_MAT_ZEROS_BASE(size.width_, size.height_, type)
#define TEST_SCALAR_MAT(){\
		Mat mat(3, 3, CV_8UCC3, Scalar(0, 0, 255));\
		std::cout << mat;\
}\

#define TEST_MAT_AT(){\
	Mat mat=Mat::zeros(3,3,CV_8SCC1);\
	mat.at<uint8_t>(2,2)=3;\
	std::cout << mat;\
}\

#define TEST_MAT_AT1(){\
	Mat mat=Mat::zeros(3,3,CV_8SCC1);\
	mat.at<uint8_t>(Point(2,2))=3;\
	std::cout << mat;\
}\


static void test_mat() {
//  Mat mat(2, 2, CV_8SC1);
//  EXPECT_EQ_INT(2, mat.cols());
//  EXPECT_EQ_INT(2, mat.rows());
// 	EXPECT_EQ_INT(1, mat.channels());
// 	EXPECT_EQ_INT(1, mat.channel_size());
// 	EXPECT_EQ_INT(CV_8SC1, mat.type());
// 	EXPECT_EQ_INT(2 * 2, mat.size());
// 	EXPECT_TRUE(mat.is_continuous());
// 	EXPECT_FALSE(mat.empty());
//	TEST_MAT_ZEROS_BASE(4, 4, CV_8UCC1);
//	TEST_MAT_ZEROS1(Size(5, 4), CV_8UCC1);
//	TEST_SCALAR_MAT();
//	TEST_MAT_AT();
//	TEST_MAT_AT1();
}

static void test() {
//	test_point();
//	test_size();
	test_mat();
}

int main() {
	test();
//	printf("%d/%d (%3.2f%%) passed\n", test_pass, test_count, test_pass * 100 / test_count);
	return 0;
}





