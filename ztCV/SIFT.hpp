#ifndef _ZT_CV_SIFT_H_
#define _ZT_CV_SIFT_H_

#include "Mat.h"
#include "Traits.h"
#include "ImageProcessing.hpp"
#include <vector>
namespace ztCV {
	
	template<typename Type>
	Mat_<Type>& initialize_image(const Mat_<Type>& src, bool scale_double, float sigma);

	template<typename Type>
	void gaussian_pyramid(const Mat_<Type>& src, std::vector<Mat_<Type>>& pyramid_vector, int octaves);




	

	template<typename Type>
	Mat_<Type>& initialize_image(const Mat_<Type>& src, bool scale_double, float sigma) {
		Mat_<Type> gray;
		if (src.channels() == 1 || src.channels() == 4) {
			rgb2gray(src, gray);
		} else {
			gray = src;
		}
		float sig_diff;
		if (scale_double) {
			sig_diff = std::sqrtf(std::pow(sigma, 2) - std::pow(SIFT_INITIAL_SIGMA, 2));
			Mat_<Type> scaled_gray(gray.rows() * 2, gray.cols() * 2, gray.type());
			resize(gray, scaled_gray, interpolation_type::INTERPOLATION_BILINEAR);
			gaussian_blur(scaled_gray, scaled_gray, Size(), sig_diff, sig_diff);
			return scaled_gray;
		} else {
			sig_diff = std::sqrtf(std::max(std::pow(sigma, 2) - std::pow(SIFT_INITIAL_SIGMA, 2)), 0.01f);
			gaussian_blur(gray, gray, Size(), sig_diff, sig_diff);
			return gray;
		}
	}

	template<typename Type>
	void gaussian_pyramid(const Mat_<Type>& src, std::vector<Mat_<Type>>& pyramid_vector, int octaves) {
		std::vector<double> sig()
	}
}

#endif