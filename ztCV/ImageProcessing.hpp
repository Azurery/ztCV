#ifndef _ZT_CV_IMAGE_PROCESSING_H_
#define _ZT_CV_IMAGE_PROCESSING_H_

#include "Mat.h"
#include "Types.h"
#include "Traits.h"
#include <cassert>
#include <functional>

namespace ztCV {

	template<typename Type>
	void box_filter(const Mat_<Type>& src, Mat_<Type>& dest, Size size, bool normalize, int border_type);

	Matf get_gaussian_kernel(Size kernel_size, double sigma_x, double sigmal_y, int kernel_type);

	template<typename Type>
	void gaussian_blur(Mat_<Type>& src, Mat_<Type>& dest, Size kernel_size, double sigma_x, double sigma_y, int border_type);

	int process_border(int pos, int width, int border_type);
	
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