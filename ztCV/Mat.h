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

/** ͼ����ȣ�
*	@concept	��ָ�洢ÿ���������õ�λ����Ҳ����������ͼ���ɫ�ʷֱ��ʵġ�
*				ȷ���˲�ɫͼ���ÿ�����ؿ����е�ɫ����������ȷ���Ҷ�ͼ���ÿ�����ؿ����еĻҶȼ�����
*				������ɫ��ͼ���п��ܳ��ֵ�����ɫ���������߻Ҷ�ͼ���е����Ҷȵȼ���
*				
*	@example	�������ص��ɫ����ϸ�ȣ���16λ(65536ɫ)��32λ�ȡ�����һ����ɫͼ����ÿ��������8λ�������Ҷ���ĿΪ256��
*				һ����ɫͼ���ÿ������������λ���ֱ�Ϊ��4��4��2���������ɫ��ĿΪ1024��
*				����˵���ص����Ϊ10λ,ÿ�����ؿ�����1024����ɫ�е�һ��.
*
*				���磺һ�����ĳߴ���1024*768�����Ϊ16��������������Ϊ1.5M��
*
*				�������£�1024*768*16bit
*						=(1024*768*16)/8�ֽ�
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
		// \parameter:	rows ��ά���������
		// \parameter:	cols ��ά���������
		// \parameter:	type ��ʾ�������͡�������ΪCV_<bit_depth>(S|U|F)C<number_of_channels>
		//						CV_8UC1 ��ָһ��8λ���޷��š����͵�ͨ���������У�ͨ����ʾÿ�����ܴ�Ŷ��ٸ�����
		//					������RGB��ɫͼ�е�ÿ�����ص�������ֵ������ͨ���ġ�
		//						8U����ʾͼƬ�е���ȱ�ʾÿ��ֵ�ɶ���λ���洢����һ���������⣬һ��ͼƬ��8bit�ģ�
		//						�������8����ʾ
		//************************************
		Mat_(int rows, int cols, int type);

		Mat_(int rows, int cols, int type, Type value);

		//*********************************** 	
		// \parameter:	size ��ά���飺Size(cols, rows)
		// \parameter:	type ��ά���������
		//			
		//************************************
		Mat_(Size size, int type);

		//************************************ 	
		// \parameter:	s ����һ����ѡֵ�����ڳ�ʼ��ÿ��Mat_Ԫ��
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
		// \brief:		���ڷ�������Ԫ�ص�����������������ʾͼ�����������	
		//************************************
		size_t total() const;
		size_t total(int dim_begin, int dim_end) const;
		size_t nums() const;

		//************************************
		// \method name:at
		// \parameter:	row ָ��Ԫ�ص���
		// \parameter:	col ָ��Ԫ�ص���
		//
		// \brief:		���ڷ���Mat������Ϊ (row,col) �ĵ�
		//				
		//				���ڶ�ά�������ڴ����԰������У��������ַΪ��
		//					addr(i,j) = data_ptr_ + layer_[0]*i + layer_[1]*j	
		//************************************
		template<typename Type2> Type2& at(int row, int col);
		template<typename Type2> const Type2& at(int row, int col) const;
		
		//************************************
		// \method name:at
		// \parameter:	row ָ��һά���������
		//
		// \brief:		������һά���飬��ֱ�ӷ���layer_[0]����	
		//************************************
		template<typename Type2> const Type2& at(int row) const;
		template<typename Type2> Type2& at(int row);

		Type& at(int row, int col);

		//************************************
		// \method name:at
		// \parameter:	pt ����Point(row,col)��ΪԪ�ص�λ�ã���ֱ������һ����
		//
		// \brief:		ֱ�Ӵ���һ���㣬������Ԫ��λ��	
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
				// FIXME: scalar�����������
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
		// \parameter:	rows ����ͼ�������
		// \parameter:	cols ����ͼ�������
		//
		// \brief:		��Mat������Ԫ�ض�����Ϊ0	
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


		// ����ָ��
		uint8_t* data_ptr_;
private:
		/** �ܹ�Ϊ13λ
		* - subMat_ flag����12λ
		* - continuity flag: ��11λ
		* - number of channels: 3-10λ
		* - depth: 0-2λ
		*/
		int flags_;
		uint32_t rows_, cols_;

		// ÿ�����ص��ͨ������
		// ���㹫ʽ��channels_ = element_size_ / channel_size_
		uint32_t channels_;
		// �����ά��
		uint8_t dims_;

		// Mat_���ݵ���ʼ��ַ��ĩβ
		uint8_t* data_begin_ = nullptr;
		uint8_t* data_end_ = nullptr;

		

		// layer_�ڴ˴�ָ����ͼ���ڸ����㼶�ϵ��ֽ����Ĵ�С
		//    �㼶��ָ����ͼ�������Ρ��������ڴ��У��Զ�άͼ��Ϊ�����ǰ��н��д洢��
		//    
		//    ĳ����(i,j)��λ��addr(i,j)���㹫ʽ��
		//			addr(i,j) = data_ptr_ + layer_[0]*i + layer_[1]*j
		//  @example ����άͼ��Ϊ��
		//	-	����һ�����棨plain����һ����layer_[0]����ʾ���ж��ٸ�plain������
		//	-	ÿһ��ƽ����һ���У�line���ڶ�����layer[1]������
		//	-	ÿ����һ�����㣨point����������layer_[2]������
		//	

		// Mat_�е� layer_[0] ����ÿһ����һ�������ڴ���ռ�ݵ��ֽ�����
		//   ���磬��άͼ���е� layer_[0] ����ÿһ�У���һ�����ھ����ڴ��У�ռ�ݵ��ֽڵ�����
		// uint32_t* layer_=nullptr;
		// 
		// layer_[0] = cols_ * element_size_
		// layer_[1] = channels_ * channel_size_ = element_size_
		// layer_(0) = layer[0] / channel_size_
		// 
		std::vector<uint32_t> layer_;


		// Mat_�е������ص���ռ���ֽ���(��ͨ�����й�)
		// ���㹫ʽ��element_size_ = channels * channel_size_
		uint8_t element_size_;

		// ���ڼ���ÿ��channel��ռ���õ��ֽ�������ͨ���޹أ�
		// ���㹫ʽ��channel_size_ = sizeof(Type)
		uint8_t channel_size_;
		
		// Mat_�ܹ�����������
		uint32_t size_;
		// ָ����һ���������������������¼�Ŷ��ٸ�Mat_ָ����ͬһ��data��������Ϊ0ʱ���ͷ�data

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
