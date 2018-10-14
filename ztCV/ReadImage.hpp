#ifndef _ZT_CV_READIMAGE_H_
#define _ZT_CV_READIMAGE_H_

// #ifdef _MSC_VER
// #define _CRT_SECURE_NO_WARNINGS
// #endif

#ifdef __cplusplus
extern "C" {
#include "jpeglib.h"
#include <setjmp.h>
}
#endif

#pragma comment(lib, "libjpeg.lib") 
#include <string>
#include "Traits.h"

namespace ztCV {
	 Mat read_image(const char file_name[]) {
		Mat mat;
		struct jpeg_decompress_struct cinfo;
		struct jpeg_error_mgr err;
		JSAMPARRAY buffer;		//输出缓存
		int row_stride;			//输出缓存的行的宽度

		// 初始化错误处理函数
		cinfo.err = jpeg_std_error(&err);
		FILE* input_file;
		if ((input_file = fopen(file_name, "rb")) == NULL) {
			fprintf(stderr, "cannot open %s\n", file_name);
			exit(1);
		}

		// 初始化解压
		jpeg_create_decompress(&cinfo);

		// 指明输入路径
		jpeg_stdio_src(&cinfo, input_file);

		// 利用jpeg_read_header()读取图片的信息，并存到mat中
		jpeg_read_header(&cinfo, TRUE);

		int type = cinfo.num_components > 1 ? CV_8UCC3 : CV_8UCC1;
		mat.create(cinfo.image_height, cinfo.image_width, type);
		//mat.dims_=

		// 开始解压
		jpeg_start_decompress(&cinfo);
		
		row_stride = cinfo.output_width*cinfo.output_components;
		buffer = (*(cinfo.mem->alloc_sarray))
			((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);
		
		while (cinfo.output_scanline < cinfo.output_height) {
			jpeg_read_scanlines(&cinfo, buffer, 1);
			// 获取到图片的像素数据：
			//		buffer保存了读取的当前行的数据，保存顺序是R、G、B
			//		output_scanline是已经读取过的行数
			int row = (cinfo.output_scanline - 1);
			for (int i = 0; i < mat.cols() ; i++) {
				mat.data_ptr_[row * mat.cols() + i] = buffer[0][i];
			}
			if (row != (mat.rows() - 2))
				row++;
		}

		jpeg_finish_decompress(&cinfo);
		jpeg_destroy_decompress(&cinfo);
		fclose(input_file);
		return mat;
	} 
}
 
#endif
	
