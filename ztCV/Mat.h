#ifndef _ZT_CV_Mat__H_
#define _ZT_CV_Mat__H_

#include <cstdint>
#include "Types.h"
#include "Traits.h"
#include <memory>
#include <vector>
#include <iterator>
#include <iostream>
#include <iomanip>

/** 图像深度：
*	@concept	是指存储每个像素所用的位数，也是用来度量图像的色彩分辨率的。
*				确定了彩色图像的每个像素可能有的色彩数，或者确定灰度图像的每个像素可能有的灰度级数。
*				决定了色彩图像中可能出现的最多的色彩数，或者灰度图像中的最大灰度等级。
*				
*	@example	单个像素点的色彩详细度，如16位(65536色)，32位等。比如一幅单色图像，若每个像素有8位，则最大灰度数目为256。
*				一幅彩色图像的每个分量的像素位数分别为：4、4、2，则最大颜色数目为1024。
*				就是说像素的深度为10位,每个像素可以是1024种颜色中的一种.
*
*				例如：一幅画的尺寸是1024*768，深度为16，则它的数据量为1.5M。
*
*				计算如下：1024*768*16bit
*						=(1024*768*16)/8字节
*
*/
namespace ztCV {

	template<typename Type>
	class Mat_ {
	public:
		void create(int rows, int cols, int type);
		
		void create(Size size, int type);
		void create(int dims, const uint32_t* size, int type);
		Mat_() = default;

		//************************************	
		// \parameter:	rows 二维数组的行数
		// \parameter:	cols 二维数组的列数
		// \parameter:	type 表示数组类型。其类型为CV_<bit_depth>(S|U|F)C<number_of_channels>
		//						CV_8UC1 是指一个8位、无符号、整型单通道矩阵。其中，通道表示每个点能存放多少个数，
		//					类似于RGB彩色图中的每个像素点有三个值，即三通道的。
		//						8U：表示图片中的深度表示每个值由多少位来存储，是一个精度问题，一般图片是8bit的，
		//						则深度是8，表示
		//************************************
		Mat_(int rows, int cols, int type);

		Mat_(int rows, int cols, int type, Type value);

		//*********************************** 	
		// \parameter:	size 二维数组：Size(cols, rows)
		// \parameter:	type 二维数组的类型
		//			
		//************************************
		Mat_(Size size, int type);

		//************************************ 	
		// \parameter:	s 它是一个可选值，用于初始化每个Mat_元素
		//		
		//************************************
		Mat_(int rows, int cols, int type, const Scalar& s);

		Mat_(Size size, int type, const Scalar& s);
		Mat_(std::initializer_list<uint8_t>& pixel_list);
		Mat_(int dims, const uint32_t* arr, int type);
		
		Mat_(int dims, const uint32_t* arr, int type, const Scalar& s);
		//Mat_(const std::vector<uint32_t>& vec, int type);
		
		Mat_(const Mat_& m, const std::vector<int>& row, const std::vector<int>& col);
		Mat_(const Mat_& m);

		~Mat_();

		Mat_& operator=(const Mat_& m);
		Mat_& operator=(const Scalar& s);
		Mat_& operator=(const std::initializer_list<int8_t>& value);

		void free();
		int add_reference_counter(int* counter, int num);
		int channels() const;
		int type() const;
		int depth() const;
		size_t element_size() const;
		size_t channel_size() const;
		bool is_continuous() const;
		bool is_submat() const;
		int rows() const;
		int cols() const;
		Mat_ row(int index) const;
		Mat_ col(int index) const;
		
		//************************************
		// \method name:total
		//
		// \brief:		用于返回数组元素的总数（如果该数组表示图像的像素数）	
		//************************************
		size_t total() const;
		size_t total(int dim_begin, int dim_end) const;
		size_t nums() const;

		//************************************
		// \method name:at
		// \parameter:	row 指定元素的行
		// \parameter:	col 指定元素的列
		//
		// \brief:		用于访问Mat中坐标为 (row,col) 的点
		//				
		//				由于二维数组在内存中以按行排列，所以其地址为：
		//					addr(i,j) = data_ptr_ + layer_[0]*i + layer_[1]*j	
		//************************************
		template<typename Type2> Type2& at(int row, int col);
		template<typename Type2> const Type2& at(int row, int col) const;
		
		//************************************
		// \method name:at
		// \parameter:	row 指定一维数组的坐标
		//
		// \brief:		由于是一维数组，则直接访问layer_[0]即可	
		//************************************
		template<typename Type2> const Type2& at(int row) const;
		template<typename Type2> Type2& at(int row);

		Type& at(int row, int col);

		//************************************
		// \method name:at
		// \parameter:	pt 利用Point(row,col)作为元素的位置，即直接输入一个点
		//
		// \brief:		直接传入一个点，返回其元素位置	
		//************************************
		template<typename Type2> const Type2& at(Point pt) const;
		template<typename Type2> Type2& at(Point pt);

 		template<int n> const Type& at(Vec_<int, n>& vec) const;
		//template<typename Type, int n> Type& at(Vec_<int, n>& vec);

		bool empty() const;

		friend std::ostream& operator<<(std::ostream& out, const Mat_<Type>& m) {
			out << "[" << std::endl;
			int size = 0;
			if (m.channels() > 1) {
				// FIXME: scalar的输出有问题
				for (int i = 0; i < m.rows(); i++) {
					for (int j = 0; j < m.cols(); j++) {
						for (int z = 0; z < 3; z++) {
							out << ' ' << *(m.data_ptr_ + z + size) << ' ';
						}
						out << ',';
						size += 3;
					}
					out << std::endl;
				}
			} else {
				for (int i = 0; i < m.rows(); i++) {
					for (int j = 0; j < m.cols(); j++) {
						out << std::setw(2) << *(m.data_ptr_ + size++) << ',';
					}
					out << std::endl;
				}
				
			}
			out << "]" << std::endl;
			return out;
		}
		
		Type* operator[](size_t n);
		const Type* operator[](size_t n) const;

		const Mat_<Type>& clone() const;

		//************************************
		// \method name:zeros
		// \parameter:	rows 输入图像的行数
		// \parameter:	cols 输入图像的列数
		//
		// \brief:		将Mat的所有元素都设置为0	
		//************************************
		static Mat_ zeros(int rows, int cols, int type);
		static Mat_ zeros(Size size, int type);

		static Mat_ ones(int rows, int cols, int type);
		static Mat_ ones(Size size, int type);

		void set_value(Type value);
		template<typename Type2> void set_value(const Scalar& s);

		void resize(Mat_<Type>& src, Mat_<Type>& dest, interpolation_type inter_type = interpolation_type::INTERPOLATION_BILINEAR);
		const Size size() const;

		Type* ptr(int row);
		const Type* ptr(int row) const;

		template<typename Type2> Type2* ptr(int row);
		template<typename Type2> const Type2* ptr(int row) const;

		Type* ptr(int row, int col);
		const Type* ptr(int row, int col) const;
		
		template<typename Type2> Type2* ptr(int row, int col);
		template<typename Type2> const Type2* ptr(int row, int col) const;

		void interpolation_nearest(Mat_<Type>& src, Mat_<Type>& dest);


		// 数据指针
		uint8_t* data_ptr_;
private:
		/** 总共为13位
		* - subMat_ flag：第12位
		* - continuity flag: 第11位
		* - number of channels: 3-10位
		* - depth: 0-2位
		*/
		int flags_;
		uint32_t rows_, cols_;

		// 每个像素点的通道数量
		// 计算公式：channels_ = element_size_ / channel_size_
		uint32_t channels_;
		// 矩阵的维数
		uint8_t dims_;

		// Mat_数据的起始地址和末尾
		uint8_t* data_begin_ = nullptr;
		uint8_t* data_end_ = nullptr;

		

		// layer_在此处指的是图像在各个层级上的字节数的大小
		//    层级是指构成图像的名层次。由于在内存中，以二维图像为例，是按行进行存储的
		//    
		//    某个点(i,j)的位置addr(i,j)计算公式：
		//			addr(i,j) = data_ptr_ + layer_[0]*i + layer_[1]*j
		//  @example 以三维图像为例
		//	-	它由一个个面（plain，第一级即layer_[0]，表示共有多少个plain）构成
		//	-	每一个平面由一行行（line，第二级即layer[1]）构成
		//	-	每行由一个个点（point，第三级即layer_[2]）构成
		//	

		// Mat_中的 layer_[0] 就是每一个第一级，在内存中占据的字节数量
		//   例如，二维图像中的 layer_[0] 就是每一行（第一级）在矩阵内存中，占据的字节的数量
		// uint32_t* layer_=nullptr;
		// 
		// layer_[0] = cols_ * element_size_
		// layer_[1] = channels_ * channel_size_ = element_size_
		// layer_(0) = layer[0] / channel_size_
		// 
		std::vector<uint32_t> layer_;


		// Mat_中单个像素点所占的字节数(与通道数有关)
		// 计算公式：element_size_ = channels * channel_size_
		uint8_t element_size_;

		// 用于计算每个channel所占用用的字节数（与通道无关）
		// 计算公式：channel_size_ = sizeof(Type)
		uint8_t channel_size_;
		
		// Mat_总共的像素数量
		uint32_t size_;
		// 指向了一个计数器，这个计数器记录着多少个Mat_指向了同一个data，当计数为0时，释放data

		int* reference_counter_ = nullptr;
		//std::shared_ptr<int32_t*> reference_counter_ = nullptr;
	};

	using Matuc = Mat_<uint8_t>;
	using Matsc = Mat_<int8_t>;
	using Matus = Mat_<uint16_t>;
	using Matss = Mat_<int16_t>;
	using Matui = Mat_<uint32_t>;
	using Matsi = Mat_<int32_t>;
	using Matf = Mat_<float>;
	using Matd = Mat_<double>;
	using Mat = Matuc;

	template<typename Type>
	Mat_<Type>& operator-(Mat_<Type>& src, Mat_<Type>& dest);


	template<typename Type>
	class mat_const_iterator {
	public:
		using value_type = Type;
		using difference_type = ptrdiff_t;
		using pointer = value_type * ;
		using reference = value_type & ;
		using iterator = pointer;
		using iterator_categoty = typename std::iterator_traits<Type>::random_access_iterator_tag;

		mat_const_iterator();
		mat_const_iterator(const Mat_<Type>& m);
		mat_const_iterator(const mat_const_iterator<Type>& other);

		const iterator begin() const;
		const iterator end() const;
		mat_const_iterator<Type>& operator=(const mat_const_iterator<Type>& other);

		const Type& operator*() const;
		const Type& operator[](difference_type diff) const;

		mat_const_iterator& operator+=(difference_type delta);
		mat_const_iterator& operator-=(difference_type delta);
		mat_const_iterator& operator--();
		mat_const_iterator operator--(int);
		mat_const_iterator& operator++();
		mat_const_iterator operator++(int);

		bool operator==(const mat_const_iterator<Type>& iter) const;
		bool operator!=(const mat_const_iterator<Type>& iter) const;
		bool operator<(const mat_const_iterator<Type>& iter) const;
		bool operator>(const mat_const_iterator<Type>& iter) const;

	private:
		iterator begin_;
		iterator end_;
		const Mat_<Type>& mat_;
		iterator cur_;
		size_t element_size_;
	};
}

#endif
