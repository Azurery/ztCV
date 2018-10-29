#ifndef _ZT_CV_RWIMAGE_H_
#define _ZT_CV_RWIMAGE_H_

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
		JSAMPARRAY buffer;		//�������
		int row_stride;			//���������еĿ��

		// ��ʼ����������
		cinfo.err = jpeg_std_error(&err);
		FILE* input_file;
		if ((input_file = fopen(file_name, "rb")) == NULL) {
			fprintf(stderr, "cannot open file %s\n", file_name);
			exit(1);
		}

		// ��ʼ����ѹ
		jpeg_create_decompress(&cinfo);

		// ָ������·��
		jpeg_stdio_src(&cinfo, input_file);

		// ����jpeg_read_header()��ȡͼƬ����Ϣ�����浽mat��
		jpeg_read_header(&cinfo, TRUE);

		int type = cinfo.num_components > 1 ? CV_8UCC3 : CV_8UCC1;
		mat.create(cinfo.image_height, cinfo.image_width, type);
		//mat.dims_=

		// ��ʼ��ѹ
		jpeg_start_decompress(&cinfo);
		
		row_stride = cinfo.output_width*cinfo.output_components;
		buffer = (*(cinfo.mem->alloc_sarray))
			((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);
		
		while (cinfo.output_scanline < cinfo.output_height) {
			jpeg_read_scanlines(&cinfo, buffer, 1);
			// ��ȡ��ͼƬ���������ݣ�
			//		buffer�����˶�ȡ�ĵ�ǰ�е����ݣ�����˳����R��G��B
			//		output_scanline���Ѿ���ȡ��������
			if (mat.channels() == 1) {
				for (int i = 0; i < mat.cols(); i++) {
					mat.at<float>(cinfo.output_scanline - 1, i) = buffer[0][i];
				}
			} else {
				for (int i = 0; i < mat.cols(); i++) {
 					mat.at<Vec3uc>(cinfo.output_scanline - 1, i)[2] = buffer[0][i * 3 + 0];
 					mat.at<Vec3uc>(cinfo.output_scanline - 1, i)[1] = buffer[0][i * 3 + 1];
 					mat.at<Vec3uc>(cinfo.output_scanline - 1, i)[0] = buffer[0][i * 3 + 2];
					//  				mat[cinfo.output_scanline - 1][i * 3 + 2] = buffer[0][i * 3 + 0];
					//  				mat[cinfo.output_scanline - 1][i * 3 + 1] = buffer[0][i * 3 + 1];
					//  				mat[cinfo.output_scanline - 1][i * 3 + 0] = buffer[0][i * 3 + 2];
				}
			}
		}

		jpeg_finish_decompress(&cinfo);
		jpeg_destroy_decompress(&cinfo);
		fclose(input_file);
		return mat;
	} 

	 bool write_image(const char* file_name, const Mat& mat) {
		 struct jpeg_compress_struct cinfo;
		 struct jpeg_error_mgr err;
		 FILE* output_file;
		 int row_stride;

		 cinfo.err = jpeg_std_error(&err);
		 jpeg_create_compress(&cinfo);

		 if ((output_file = fopen(file_name, "wb")) == NULL) {
			 fprintf(stderr, "cannot fopen file %s\n", file_name);
			 exit(1);
		 }

		 jpeg_stdio_dest(&cinfo, output_file);
		 cinfo.image_width = static_cast<JDIMENSION>(mat.cols());
		 cinfo.image_height = static_cast<JDIMENSION>(mat.rows());
		 cinfo.input_components = mat.channels();
		 cinfo.in_color_space = (mat.channels() == 1 ? JCS_GRAYSCALE : JCS_RGB);

		 //�趨ȱʡ����
		 jpeg_set_defaults(&cinfo);

		 //quality�Ǹ�0��100֮�����������ʾѹ�����ʡ�
		 int quality = 40;
		 jpeg_set_linear_quality(&cinfo, quality, TRUE);

		 jpeg_start_compress(&cinfo, TRUE);

		 // д������
		 if (mat.channels() == 1) {
			 row_stride = mat.cols();
		 } else {
			 row_stride = mat.cols() * 3;
		 }
		 JSAMPROW row_pointer[1];

		 // ѹ������ʹ�õ�cinfo.next_scanline < cinfo.image_height���ж��Ƿ����д�����ݡ�
		 while (cinfo.next_scanline < cinfo.image_height) {
			 //�ҵ�ͼ���е�ĳһ�У�д��Ŀ���ļ�
//   			 uint8_t* tmp = (uint8_t*)malloc(sizeof(uint8_t)*row_stride);
//   			 for (int i = 0; i < row_stride; i++) {
//   				tmp[i] = mat.data_ptr_[cinfo.next_scanline * row_stride + i];
//   			 }
			 row_pointer[0] = &mat.data_ptr_[cinfo.next_scanline * row_stride];
			 (void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
		 }

		 jpeg_finish_compress(&cinfo);
		 fclose(output_file);
		 jpeg_destroy_compress(&cinfo);
	 }
}
 
#endif
	
