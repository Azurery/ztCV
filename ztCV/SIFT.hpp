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

		Mat_<Type>& derivation_3d(std::vector<std::vector<Mat_<Type>>>& dog_pyramid, int octave, int interval, int r, int c);

		bool adjust_local_extremum(const std::vector<std::vector<Mat_<Type>>>& dog_pyramid, int octave, int interval, int r, int c, int intervals);
		
		double calculate_orientation_histogram(const Mat_<Type>& img, Point point, int radius, double sigma);
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
	
		for (int o = 0; o < octaves; o++) {
			for (int i = 1; i <= intervals; i++) {
				const Mat_<Type>& cur = dog_pyramid[o][i];
				const Mat_<Type>& pre = dog_pyramid[o][i - 1];
				const Mat_<Type>& next = dog_pyramid[o][i + 1];
				int rows = cur.rows(), cols = cur.cols();
				int layer = cur.layer_[0];
				// ���ӣ�5+3���G�����5+2��DOG�㣬Ȼ��ȥͷȥβ�ӵ�1����6�㣬��6��ı�ž���octaves
				// SIFT_IMG_BORDER = 5; �������ֻ��Ϊ�˷�ֹԽ��
				// ��5��ʼ��rows-5����ֹԽ��
				for (int r = SIFT_IMAGE_BORDER; r < rows - SIFT_IMAGE_BORDER; r++) {
					auto cur_ptr = cur.ptr(r);
					auto pre_ptr = pre.ptr(r);
					auto next_ptr = next.ptr(r);
					for (int c = SIFT_IMAGE_BORDER; c < cols - SIFT_IMAGE_BORDER; c++) {
						// �ж��Ƿ�Ϊ��ֵ��
						// ԭ��Ϊ:ͨ���͸�˹����������һ���9������ + ����ĳ��˱������Լ���������8������
						// ����һ���9�����ؽ��бȽϿ��Ƿ�Ϊ��26����������С��һ�������Ƿ�Ϊ����һ����
						// �������Ϊ��ֵ�㡣
						int value = cur_ptr[c];
						// //����Сֵ����Ҫ�������ϲ��²����Χ3*3���������������жϼ�ֵ
						if (std::abs(value) > threshold &&
							((	value > 0						 &&
								value >= cur_ptr[c - 1]			 && value >= cur_ptr[c + 1]		 &&
								value >= cur_ptr[c - layer - 1]  && value >= cur_ptr[c - layer]  && value >= cur_ptr[c - layer + 1]  &&
								value >= cur_ptr[c + layer + 1]  && value >= cur_ptr[c + layer]  && value >= cur_ptr[c + layer + 1]  &&
								value >= next_ptr[c - 1]		 && value >= next_ptr[c]		 && value >= next_ptr[c + 1]		 &&
								value >= next_ptr[c - layer - 1] && value >= next_ptr[c - layer] && value >= next_ptr[c - layer + 1] &&
								value >= next_ptr[c + layer - 1] && value >= next_ptr[c + layer] && value >= next_ptr[c + layer + 1] &&
								value >= pre_ptr[c - 1]			 && value >= pre_ptr[c]			 && value >= pre_ptr[c + 1]			 &&
								value >= pre_ptr[c - layer - 1]  && value >= pre_ptr[c - layer]  && value >= pre_ptr[c - layer + 1]  &&
								value >= pre_ptr[c + layer - 1]  && value >= pre_ptr[c + layer]  && value >= pre_ptr[c + layer + 1]) ||
							(	value < 0						 &&
								value <= cur_ptr[c - 1]			 && value <= cur_ptr[c + 1]		 &&
								value <= cur_ptr[c - layer - 1]  && value <= cur_ptr[c - layer]  && value <= cur_ptr[c - layer + 1]  &&
								value <= cur_ptr[c + layer + 1]  && value <= cur_ptr[c + layer]  && value <= cur_ptr[c + layer + 1]  &&
								value <= next_ptr[c - 1]		 && value <= next_ptr[c]		 && value <= next_ptr[c + 1]		 &&
								value <= next_ptr[c - layer - 1] && value <= next_ptr[c - layer] && value <= next_ptr[c - layer + 1] &&
								value <= next_ptr[c + layer - 1] && value <= next_ptr[c + layer] && value <= next_ptr[c + layer + 1] &&
								value <= pre_ptr[c - 1]			 && value <= pre_ptr[c]			 && value <= pre_ptr[c + 1]			 &&
								value <= pre_ptr[c - layer - 1]  && value <= pre_ptr[c - layer]  && value <= pre_ptr[c - layer + 1]  &&
								value <= pre_ptr[c + layer - 1]  && value <= pre_ptr[c + layer]  && value <= pre_ptr[c + layer + 1]))) {
							int r1 = r, c1 = c, interval = i;

						}
					}
				}
			}
		}
	}

	//������άƫ����
	// ������x��y�����ϵ�ƫ��������˹��ֳ߶ȿռ�����������صĳ߶�
	// ʵ��������ɢ�����м���ƫ������ͨ���������ص�����������
	// ����˵����x�����ƫ����dx����ͨ����������x����ĺ�һ����ȥǰһ����Ȼ�����2���������dx
	template<typename Type>
	Mat_<Type>& ztCV::SIFT<Type>::derivation_3d(std::vector<std::vector<Mat_<Type>>>& dog_pyramid, int octave, int interval, int r, int c) {
		double dx, dy, ds;
		const Mat_<Type>& cur = dog_pyramid[octave][interval];
		const Mat_<Type>& pre = dog_pyramid[octave][interval - 1];
		const Mat_<Type>& next = dog_pyramid[octave][interval + 1];

		dx = (cur[r][c + 1] - cur[r][c - 1]) / 2.0;
		dy = (cur[r + 1][c] - cur[r - 1][c]) / 2.0;
		ds = (next[r][c] - pre[r][c]) / 2.0;
		Matd dirivation(Size(3, 1), CV_64DC1);
		dirivation = { dx,dy,ds };
		return dirivation;
	}

	template<typename Type>
	bool SIFT<Type>::adjust_local_extremum(const std::vector<std::vector<Mat_<Type>>>& dog_pyramid, int octave, int interval, int r, int c, int intervals) {
		
		int i = 0;
		double xi = 0, xr = 0, xc = 0;
		while( i < SIFT_MAX_INTERPOLATION_STEPS) {
			const Mat_<Type>& cur = dog_pyramid[octave][interval];
			const Mat_<Type>& pre = dog_pyramid[octave][interval - 1];
			const Mat_<Type>& next = dog_pyramid[octave][interval + 1];
			double cur_value = static_cast<double>(cur[r][c]);
			Matd derivation = derivation_3d(dog_pyramid, octave, interval, r, c);
			// ������ά��ɭ������
			double dxx = (cur[r][c + 1] + cur[r][c - 1] - 2 * cur_value);
			double dyy = (cur[r + 1][c] + cur[r - 1][c] - 2 * cur_value);
			double dss = (next[r + 1][c+1] + pre[r - 1][c] - 2 * cur_value);
			double dxy = (cur[r + 1][c + 1] + cur[r + 1][c - 1] -
						  cur[r - 1][c + 1] + cur[r - 1][c - 1]) / 4.0;
			double dxs = (cur[r][c + 1] - cur[r][c - 1] -
						  pre[r, c + 1] + pre[r - 1][c]) / 4.0;
			double dys = (next[r + 1][c] - next[r - 1][c] -
						  pre[r + 1][c] + pre[r - 1][c]) / 4.0;

			Matf hessian(3, 3, CV_64DC1);
			hessian = { dxx, dxy, dxs,
						dxy, dyy, dys,
						dxs, dys, dss };
			// FIXME
			Vec3d x = hessian.solve(derivation);
			xi = -x[2];
			xr = -x[1];
			xc = -x[0];
			
			if (std::abs(xi) < 0.5&&std::abs(xr) < 0.5&&std::abs(xc) < 0.5) {
				break;
			}
			c += std::round(xc);
			r += std::round(xr);
			interval += std::round(xi);
			
			if (interval<1 ||
				interval>intervals ||
				c < SIFT_IMAGE_BORDER ||
				r < SIFT_IMAGE_BORDER ||
				c >= cur.cols() - SIFT_IMAGE_BORDER ||
				r >= cur.rows() - SIFT_IMAGE_BORDER)
				return false;

			i++;
		}

		//ȷ����ֵ���Ǿ������5���ҵ���
		if (i >= SIFT_MAX_INTERPOLATION_STEPS) {
			return false;
		}
		
		// ��ȡ�ҵ��ļ�ֵ��ĶԱȶ�
		Matd differitial = derivation_3d(dog_pyramid, octave, interval, r, c);

	}


	template<typename Type>
	double SIFT<Type>::calculate_orientation_histogram(const Mat_<Type>& img, Point point, int radius, double sigma) {
		int i, j, length = std::pow(radius * 2 + 1, 2);
		//double magnitude, oritention;
		std::vector<double> magnitude, oritention;
		// ���������ؼ���Ҷȷ���ֱ��ͼ
		// �Թؼ���Ϊ���ĵ����򴰿��ڲ���������ֱ��ͼͳ���������ص��ݶ�
		// �����ݶ�ֱ��ͼ�ķ�Χ��0��360�ȣ�����ÿ10��һ�������ܹ�36����
		for (i = -radius; i <= radius; i++) {
			for (j = -radius; j <= radius; j++) {
				int x = point.x_ + i;
				int y = point.y_ + j;
				if (x > 0 && x < img.cols() && y>0 && y < img.rows()) {
					double dx = img[x][y + 1] - img[x][y - 1];
					double dy = img[x - 1][y] - img[x + 1][y];
					magnitude.push_back(std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)));
					oritention.push_back(std::atan2(dy, dx));
				}
			}
		}
	}


}

#endif