#ifndef _ZT_CV_IMAGE_PROCESSING_H_
#define _ZT_CV_IMAGE_PROCESSING_H_

#include "Mat.h"
#include "Types.h"
#include "Traits.h"
#include <cassert>

namespace ztCV {

	template<typename Type>
	void box_filter(const Mat_<Type>& src, Mat_<Type>& dest, int dest_depth, Size size, Point point, bool normalize, int border_type);

	Matf get_gaussian_kernel(Size kernel_size, double sigma_x, double sigmal_y, int kernel_type);
	
	template<typename Type>
	void gaussian_blur(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size, double sigma_x, double sigma_y, int border_type);

	int process_border(int pos, int width, int border_type);


	///////////////// implementation ///////////////

	//////////////// box_filter ///////////////////
	template<typename Type>
	void box_filter(const Mat& src, Mat& dest, int dest_depth, Size size, Point point, bool normalize, int border_type) {
		int src_type = src.type();
		int src_depth = src.depth();
		int src_channel = src.channels();
		if (dest_depth < 0) {
			dest_depth = src_depth;
		}
		int dest_type = MAKE_TYPE(dest_depth, src_channel);
		dest.create(dest.rows(), dest.cols(), dest_type);
		//if(border_type)
	}

	Matf get_gaussian_kernel(Size kernel_size, double sigma_x, double sigma_y, int kernel_type) {
// 		assert(kernel_type == static_cast<int>(cv_depth::CV_FLOAT)
// 			|| kernel_type == static_cast<int>(cv_depth::CV_DOUBLE));
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
				gaussian_kernel.at<Vec3f>(i, j) = tmp;
				sum += tmp;
			}
		}

//    		for (int i = 0; i < gaussian_kernel.rows(); i++) {
//    			for (int j = 0; j < gaussian_kernel.cols(); j++) {
//    				std::cout << gaussian_kernel.at<float>(i, j) << " ";
//    				
//    				std::cout << ",";
//    			}
//    			std::cout << std::endl;
//    		}

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
		auto channels = src.channels();
		dest.create(src.rows(), src.cols(), src.type());
		auto tmp = new float[channels];
		
		Matf gaussian_kernel = get_gaussian_kernel(kernel_size, sigma_x, sigma_y, CV_32FC3);
		for (int i = 0; i < dest.rows(); i++) {
			for (int j = 0; j < dest.cols(); j++) {
				memset(tmp, 0, channels);
				
				for (int k = 0; k < kernel_size.width_; k++) {
					for (int l = 0; l < kernel_size.height_; l++) {
						
						auto new_i = i - kernel_size.width_ / 2 + k;
						auto new_j = j - kernel_size.height_ / 2 + l;

						if (!(new_i < dest.rows()) || !(new_j < dest.cols()) || (new_i < 0) || (new_j < 0)) {
							new_i = process_border(new_i, dest.rows(), border_type);
							new_j = process_border(new_j, dest.cols(), border_type);
						}

						for (int m = 0; m < channels; m++) {
							tmp[k] += src.at<Vec3f>(new_i, new_j)[m] * gaussian_kernel.at<float>(k, l);
						}
					}
				}

				for (int n = 0; n < channels; n++) {
					dest.at<Vec3f>(i, j)[n] = tmp[n];
				}

			}
			std::cout << i;
		}
	
		delete[] tmp;
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
}

#endif