#ifndef _ZT_CV_IMAGE_PROCESSING_H_
#define _ZT_CV_IMAGE_PROCESSING_H_

#include "Mat.h"
#include "Types.h"
#include "Traits.h"
#include <cassert>
#include <functional>
#include <algorithm>

namespace ztCV {

	template<typename Type>
	void box_filter(const Mat_<Type>& src, Mat_<Type>& dest, Size size, bool normalize, int border_type);

	template<typename Type>
	void normalized_filter(const Mat_<Type>& src, Mat_<Type>& dest, Size size, int border_type) {
		box_filter(src, dest, size, border_type);
	}

	Matf get_gaussian_kernel(Size kernel_size, double sigma_x, double sigmal_y, int kernel_type);

	template<typename Type>
	void gaussian_blur(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size, double sigma_x, double sigma_y, int border_type);

	int process_border(int pos, int width, int border_type);
	
	template<typename Type>
	void median_filter(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);

	template<typename Type, typename Type2>
	void apply(const Mat_<Type>& src, Mat_<Type>& dest, Mat_<Type2>& kernel, int border_type) {
		
		int channels = src.channels();
		int border = kernel.rows() / 2;
		int rows = dest.rows() - border;
		int cols = dest.cols() - border;
		for (int i = border; i < rows; i++) {
			for (int j = border; j < cols; j++) {
				
				float sum[3] = {};
				for (int k = -border; k <= border; k++) {
					for (int l = -border; l <= border; l++) {

// 						int new_i = i - kernel.rows()/ 2 + k;
// 						int new_j = j - kernel.cols() / 2 + l;
// 
//  						if (!(new_i < dest.rows()) || !(new_j < dest.cols()) || (new_i < 0) || (new_j < 0)) {
//  							new_i = process_border(new_i, dest.rows(), border_type);
//  							new_j = process_border(new_j, dest.cols(), border_type);
//  						}
						if (channels == 1) {
							sum[0] += src.at<float>(i + k, j + l) * kernel.at<float>(border + j, border + k);
						} else {
							auto elem = kernel.at<float>(border + k, border + l);
							sum[0] += elem * src.at<Vec3uc>(i + k, j + l)[0];
							sum[1] += elem * src.at<Vec3uc>(i + k, j + l)[1];
							sum[2] += elem * src.at<Vec3uc>(i + k, j + l)[2];
						}
					}
				}
				
				for (int k = 0; k < channels; k++) {
					if (sum[k] < 0)
						sum[k] = 0;
					else if (sum[k] > 255)
						sum[k] = 255;
				}

				if (channels == 1)
					dest.at<uint8_t>(i, j) = static_cast<uint8_t>(sum[0]);
				else if (channels == 3) {
					Vec3uc rgb = { 
						static_cast<uint8_t>(sum[2]), 
						static_cast<uint8_t>(sum[1]), 
						static_cast<uint8_t>(sum[0]) 
					};
					dest.at<Vec3uc>(i, j) = rgb;
				}
			}
			std::cout << i;
		}
	}

	float rgb2gray(Vec3uc& rgb);


	///////////////// 形态学滤波 //////////////////////////
	Mat get_structuring_element(morphology_shape shape, Size kernel_size, Point anchor);

	template<typename Type>
	void morphology_operation(Mat_<Type>& src, Mat_<Type>& dest, Size size, morphology_type mor_type);

	//************************************
	// \method name:erode
	//
	// \brief:	dst(x, y) = min{ src(x + x′, y + y′),(x′, y′) : element(x′, y′)≠0 }
	// 
	//			腐蚀是求局部最小值的操作，腐蚀操作会使图像中的高亮区逐渐减小
	//************************************
	template<typename Type>
	void erode(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);
	
	
	//************************************
	// \method name:dilate
	//
	// \brief:	dst(x, y) = max{src(x + x′, y + y′),(x′, y′) :element(x′, y′)≠0}
	//			
	//			膨胀就是求局部最大值的操作。从数学角度来说，就是将图像与核进行卷积，计算核B覆盖区域的像素点的最大值，
	//			并把这个最大值赋值给参考点指定的元素。这样就会使图像中的高亮区域逐渐增长。
	//************************************
	template<typename Type>
	void dilate(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);



	//************************************
	// \method name:open
	//
	// \brief:	dst=open(src,element)=dilate(erode(src,element))
	//			先腐蚀后膨胀
	//			
	//			用于消除小物体，在纤细点处分离物体，并且在平滑较大物体的边界的同时不明显改变其面积，
	//			同时抑制比结构元小的亮细节。 
	//************************************
	template<typename Type>
	void open(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);


	//************************************
	// \method name:close
	//
	// \brief:	dst=close(src,element)=erode(dilate(src,element))
	//			先膨胀后腐蚀 
	//			
	//			用来填充物体内细小空洞、连接邻近物体、平滑其边界的同时并不明显改变其面积，
	//			同时抑制比结构元小的暗细节
	//************************************
	template<typename Type>
	void close(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);

	//************************************
	// \method name:gradient
	//
	// \brief:	膨胀图和腐蚀图相减
	// 
	//			对二值化图像进行这一操作可以将边缘突出来，可以使用形态学梯度来保留物体的边缘轮廓
	//************************************
	template<typename Type>
	void gradient(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);

	template<typename Type>
	void tophat(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);

	template<typename Type>
	void blackhat(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);

	template<typename Type>
	void hitmiss(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size);


	
	///////////////// implementation ///////////////

	//////////////// box_filter ///////////////////

	template<typename Type>
	void box_filter(const Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size, bool normalize, int border_type) {
		assert(kernel_size.width_ % 2 == 1);
		int dest_type = src.type();
		dest.create(src.rows(), src.cols(), dest_type);
		float coefficient = 1.0;
		if (normalize)
			coefficient = 1.0 / (kernel_size.width_*kernel_size.height_);
		Matf kernel(kernel_size.width_, kernel_size.height_, CV_32FC1, coefficient);
		apply(src, dest, kernel, border_type);
	}


	Matf get_gaussian_kernel(Size kernel_size, double sigma_x, double sigma_y, int kernel_type) {
 		assert(kernel_type == static_cast<int>(cv_depth::CV_FLOAT)
 			|| kernel_type == static_cast<int>(cv_depth::CV_DOUBLE));
		Matf gaussian_kernel(kernel_size, kernel_type);
// 		float* kernel_float = reinterpret_cast<float*>(gaussian_kernel.data_ptr_);
// 		double* kernel_double = reinterpret_cast<double*>(gaussian_kernel.data_ptr_);

		sigma_x = (sigma_x > 0 ? sigma_x : 0.3*((kernel_size.width_ - 1)*0.5 - 1) + 0.8);
		sigma_y = (sigma_y > 0 ? sigma_y : 0.3*((kernel_size.width_ - 1)*0.5 - 1) + 0.8);
		double scale_2x = -1 / (2 * std::pow(sigma_x, 2));
		double scale_2y = -1 / (2 * std::pow(sigma_y, 2));
		double sum = 0;

		double x_center = kernel_size.width_ / 2;
		double y_center = kernel_size.height_ / 2;
		for (int i = 0; i < gaussian_kernel.rows(); i++) {
			for (int j = 0; j < gaussian_kernel.cols(); j++) {
				double tmp = std::exp(-(scale_2x*std::pow(i - x_center, 2))) + std::exp(-(scale_2y*std::pow(j - y_center, 2)));
				gaussian_kernel.at<float>(i, j) = tmp;
				sum += tmp;
			}
		}

 		for (int i = 0; i < gaussian_kernel.rows(); i++) {
 			for (int j = 0; j < gaussian_kernel.cols(); j++) {
 				gaussian_kernel.at<float>(i, j) /= sum;
 			}
 		}

		return gaussian_kernel;
	}

	template<typename Type>
	void gaussian_blur(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size, double sigma_x, double sigma_y, int border_type) {
		assert(kernel_size.width_ % 2 == 1 && kernel_size.height_ % 2 == 1);
		dest.create(src.rows(), src.cols(), src.type());
		
		Matf gaussian_kernel = get_gaussian_kernel(kernel_size, sigma_x, sigma_y, CV_32FC1);
		apply(src, dest, gaussian_kernel, border_type);
	}

#define SORT_COLOR(color) {\
std::sort(color.begin(), color.end(), [](int a, int b) {return a < b; });\
}

 #define PROCESS_COLOR(buffer, src, red_sum, green_sum, blue_sum, i, j, k ,l, c)\
 	do{\
 		buffer.at<Vec3uc>(k + border, l + border)[c] = src.at<Vec3uc>(i + k, j + l)[c];\
 		auto pos = buffer.at<Vec3uc>(k + border, l + border);\
 		red_sum.push_back(pos[2]);\
 		green_sum.push_back(pos[1]);\
 		blue_sum.push_back(pos[0]); \
 	} while(0)

	// FIXME:
	template<typename Type>
	void median_filter(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size) {
		int border = kernel_size.width_ / 2;
		int rows = src.rows() - border;
		int cols = src.cols() - border;
		dest.create(src.rows(), src.cols(), src.type());
		Mat buffer(kernel_size.width_, kernel_size.height_, CV_8UCC3);

		for (int i = border; i < rows; i++) {
			for (int j = border; j < cols; j++) {

				std::vector<uint8_t> red_sum, green_sum, blue_sum;
				for (int k = -border; k < border; k++) {
					for (int l = -border; l < border; l++) {
						for (int c = 0; c < src.channels(); c++) {
							PROCESS_COLOR(buffer, src, red_sum, green_sum, blue_sum, i, j, k, l, c);
						}
					}
				}

				SORT_COLOR(red_sum);
				SORT_COLOR(green_sum);
				SORT_COLOR(blue_sum);

				int center = border + 1;
				Vec3uc rgb = {
					static_cast<uint8_t>(red_sum[center]),
					static_cast<uint8_t>(green_sum[center]),
					static_cast<uint8_t>(blue_sum[center])
				};
				dest.at<Vec3uc>(i, j) = rgb;
			}
			std::cout << i;
		}
	}


	int process_border(int pos, int width, int border_type) {
		switch (border_type) {
			case static_cast<int>(border_type::BORDER_CONSTANT) :
				break;
			case static_cast<int>(border_type::BORDER_REPLICATE) :
				pos = (pos < 0 ? 0 : width - 1);
				break;
			case static_cast<int>(border_type::BORDER_WRAP) :
			{
				if (pos < 0) {
					pos -= ((pos - width + 1) / width)*width;
				}
				if (pos >= width) {
					width %= width;
				}
			}
		}
		return pos;
	}

	float rgb2gray(Vec3uc& rgb) {
		return (0.29900 * rgb[0] + 0.58700 * rgb[1] + 0.11400 * rgb[2]);
	}

	Mat get_structuring_element(morphology_shape shape, Size kernel_size, Point anchor) {
		// 只支持三种形状的结构元（structuring element）
		assert(shape == morphology_shape::MORPHOLOGY_RECTANGLE ||
			shape == morphology_shape::MORPHOLOGY_CROSS ||
			shape == morphology_shape::MORPHOLOGY_ELLIPSE);

		int semimajor, minor;
		// 当kernel为一个点时，就直接看成矩形
		if (kernel_size == Size(1, 1)) {
			shape = morphology_shape::MORPHOLOGY_RECTANGLE;
		}

		if (shape == morphology_shape::MORPHOLOGY_ELLIPSE) {
			semimajor = kernel_size.width_ / 2;
			minor = kernel_size.height_ / 2;
			
		}
		// FIXEME:先忽略cross和ellipse
		Mat mat(kernel_size, CV_8UCC1);
		for (int i = 0; i < kernel_size.height_; i++) {
			int j1 = 0, j2 = 0;
			if (shape == morphology_shape::MORPHOLOGY_RECTANGLE) {
				j2 = kernel_size.width_;
				for (int j = 0; j < j2; j++) {
					mat[i][j] = 1;
				}
			}
// 			else if (shape == morphology_shape::MORPHOLOGY_CROSS) {
// 
// 			} else {
// 				
// 			}
		}
		return mat;
	}

	template<typename Type>
	void morphology_operation(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size, morphology_type mor_type) {
		dest.create(src.rows(), src.cols(), src.type());
		Mat_<Type> buffer(kernel_size, src.type());

		int border = kernel_size.width_ / 2;
		int rows = src.rows() - border;
		int cols = src.cols() - border;

		for (int i = border; i < rows; i++) {
			for (int j = border; j < cols; j++) {
				std::vector<uint8_t> red_sum, green_sum, blue_sum;
				for (int k = -border; k < border; k++) {
					for (int l = -border; l < border; l++) {
						for (int c = 0; c < src.channels(); c++) {
							PROCESS_COLOR(buffer, src, red_sum, green_sum, blue_sum, i, j, k, l, c);
						}
					}
				}

				if (mor_type == morphology_type::MORPHOLOGY_ERODE) {
					Vec3uc rgb = {
						static_cast<uint8_t>(*std::min_element(red_sum.begin(),red_sum.end())),
						static_cast<uint8_t>(*std::min_element(green_sum.begin(),green_sum.end())),
						static_cast<uint8_t>(*std::min_element(blue_sum.begin(),blue_sum.end())),
					};
					dest.at<Vec3uc>(i, j) = rgb;
				} else {
					Vec3uc rgb = {
						static_cast<uint8_t>(*std::max_element(red_sum.begin(),red_sum.end())),
						static_cast<uint8_t>(*std::max_element(green_sum.begin(),green_sum.end())),
						static_cast<uint8_t>(*std::max_element(blue_sum.begin(),blue_sum.end())),
					};
					dest.at<Vec3uc>(i, j) = rgb;
				}
			}
			std::cout << i;
		}
	}


	template<typename Type>
	void erode(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size) {
		morphology_operation(src, dest, kernel_size, morphology_type::MORPHOLOGY_ERODE);
	}

	template<typename Type>
	void dilate(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size) {
		morphology_operation(src, dest, kernel_size, morphology_type::MORPHOLOGY_DILATE);
	}

	template<typename Type>
	void open(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size) {
		dest.create(src.rows(), src.cols(), src.type());
		erode(src, dest, kernel_size);
		dilate(src, dest, kernel_size);
	}

	template<typename Type>
	void close(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size) {
		dest.create(src.rows(), src.cols(), src.type());
		dilate(src, dest, kernel_size);
		erode(src, dest, kernel_size);
	}

	template<typename Type>
	void gradient(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size) {
		auto src_img = src;
		Mat_<Type> dilate_ret, erode_ret;
		dilate(src, dilate_ret, kernel_size);
		erode(src, erode_ret, kernel_size);
		dest = dilate_ret - erode_ret;
	}



}
#endif