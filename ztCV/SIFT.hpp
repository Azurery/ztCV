#ifndef _ZT_CV_SIFT_H_
#define _ZT_CV_SIFT_H_

#include "Mat.h"
#include "Traits.h"
#include "ImageProcessing.hpp"
#include <vector>
namespace ztCV {
	
	template<typename Type>
	class SIFT {
		Mat_<Type>& initialize_image(const Mat_<Type>& src, bool scale_double, float sigma);

		void gaussian_pyramid(const Mat_<Type>& src, std::vector<std::vector<Mat_<Type>>>& gaussian_pyramid, int octaves);

		void dog_pyramid(const std::vector<std::vector<Mat_<Type>>>& gaussian_pyramid, std::vector<std::vector<Mat_<Type>>>& dog_pyramid, int octaves, int intervals);

		void find_scale_space_extremum(const std::vector<std::vector<Mat_<Type>>>& gaussian_pyramid, const std::vector<std::vector<Mat_<Type>>>& dog_pyramid,
			int octaves, int intervals, double contrast_threshold, int edge_feature_threshold);

	private:
		int intervals;
	};

	template<typename Type>
	Mat_<Type>& SIFT<Type>::initialize_image(const Mat_<Type>& src, bool scale_double, float sigma) {
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
	void SIFT<Type>::gaussian_pyramid(const Mat_<Type>& src, std::vector<std::vector<Mat_<Type>>>& gaussian_pyramid, int octaves) {
		// Ϊ����ÿ���м��S���߶ȵļ�ֵ�㣬��DoG������ÿ����S+2��ͼ��
		// DoG�������ɸ�˹������������������õ������˹������ÿ����S+3��ͼ��
		std::vector<double> sigma(intervals + 3);
		sigma[0] = sigma;
		// k == 2^(1/s)
		// \sigma_{total}^2 = \sigma_{i}^2 + \sigma_{i-1}^2
		double k = std::pow(2.f, 1.f / intervals);
		for (int i = 1; i < intervals + 3; i++) {
			// ÿһ���sigma
			double sigma_pre = std::pow(k, static_cast<double>(i - 1)) * sigma;
			double sigma_total = sigma_pre * k;
			sigma[i] = std::sqrt(std::pow(sigma_total, 2) - std::pow(sigma_pre, 2));
		}

		for (int o = 0; o < octaves; i++) {
			gaussian_pyramid.push_back(std::vector<Mat_<Type>>(intervals + 3));
		}

		for (int o = 0; o < octaves; o++) {
			for (int i = 0; i < intervals; i++) {
				Mat_<Type>& dest = gaussian_pyramid[o][i];
				// ����ǵ�0���������ĵ�0�㣬ֱ�ӽ��и�ֵ
				if (o == 0 && i == 0) {
					dest = src;
				} else if (i == 0) {
				// �������һ���������ĵ�0�㣬Ҫʹ�õ�intervals���������
				// ÿһ��ĵ�һ�㶼��ͨ����ǰ��һ���������һ��Ľ�����ʵ�ֵ�
					const Mat_<Type>& src = gaussian_pyramid[o][i];
					resize(src, dest(src.rows() / 2, src.cols()), interpolation_type::INTERPOLATION_BILINEAR);
				} else {
					// ÿһ�����������ʹͨ��ʹ�ò�ͬsigma�ĸ�˹ģ�������д���
					const Mat_<Type>& src = gaussian_pyramid[o][i - 1];
					gaussian_blur(src, dest, Size(), sigma[i], sigma[i]);
				}
			}
		}
	}


	template<typename Type>
	void SIFT<Type>::dog_pyramid(const std::vector<std::vector<Mat_<Type>>>& gaussian_pyramid, std::vector<std::vector<Mat_<Type>>>& dog_pyramid, int octaves, int intervals) {
		// DoG�������Ĳ����Ǹ�˹����������-1
		for (int i = 0; i < octaves; i++) {
			dog_pyramid.push_back(std::vector<Mat_<Type>>(intervals + 2));
		}
		
		//dog_pyramid.resize(intervals + 2);
		for (int o = 0; o < octaves; o++) {
			for (int i = 0; i < intervals + 2; i++) {
				// ȡ��ǰ��
				const Mat_<Type>& src1 = gaussian_pyramid[o][i];
				// ȡ����һ��
				const Mat_<Type>& src2 = gaussian_pyramid[o][i + 1];
				Mat_<Type>& dest = dog_pyramid[o][i];
				dest = subtract(src1 - src2);
			}
		}
	}

	// constrast_threshold��ȥ���Աȶȵ͵ĵ������õ���ֵ
	// edge_feature_threshold��ȥ����Ե��������ֵ
	template<typename Type>
	void SIFT<Type>::find_scale_space_extremum(const std::vector<std::vector<Mat_<Type>>>& gaussian_pyramid, const std::vector<std::vector<Mat_<Type>>>& dog_pyramid,
		int octaves, int intervals, double contrast_threshold, int edge_feature_threshold) {
		int threshold = std::floor(0.5 * contrast_threshold * intervals);
	}


}

#endif