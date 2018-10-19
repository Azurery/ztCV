#ifndef _ZT_CV_IMAGE_PROCESSING_H_
#define _ZT_CV_IMAGE_PROCESSING_H_

#include "Mat.h"
#include "Mat.hpp"
#include "Types.h"
#include "Traits.h"

namespace ztCV {

	template<typename Type>
	void box_filter(const Mat& src, Mat& dest, int dest_depth, Size size, Point point, bool normalize, int border_type);



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

	}
}

#endif