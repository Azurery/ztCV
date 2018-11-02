#include "Traits.h"
#include <vector>
#include <cassert>
namespace ztCV {
	template<typename Type>
	void Mat_<Type>::create(int rows, int cols, int type) {
		//assert(rows_*cols_ != 0 && "the image is valid!");

		this->rows_ = rows;
		this->cols_ = cols;
		this->channel_size_ = sizeof(Type);
		this->channels_ = (type >> 3) & 0xFF;
		this->size_ = this->rows_ * this->cols_*this->channels_;

		this->flags_ = (MAT_FLAGS_MASK & type) | (DataType<Type>::depth << 8) | (DataType<Type>::channels);
		this->element_size_ = this->channel_size_*this->channels_;
		this->layer_ = { (this->cols_*this->element_size_),this->element_size_ };


		// 如果reference counter为1，则释放内存
		//free();

		// 重新分配内存
		this->data_ptr_=this->data_begin_= 
			reinterpret_cast<uint8_t*>(new uint8_t[this->size_]);
		this->data_end_ = this->data_ptr_ + this->size_*this->element_size_;
		
		this->reference_counter_ = new int(1);
	}


	template<typename Type>
	void Mat_<Type>::create(Size size, int type) {
		create(size.width_, size.height_, type);
	}

	template<typename Type>
	int Mat_<Type>::add_reference_counter(int* counter, int num) {
		const int tmp = *counter;
		*counter += num;
		return tmp;
	}


	template<typename Type>
	void Mat_<Type>::free() {
		if (reference_counter_&&add_reference_counter(reference_counter_, -1) == 1) {
			delete[] reinterpret_cast<Type*>(data_ptr_);
			data_ptr_ = data_begin_ = data_end_ = nullptr;

			delete reference_counter_;
			reference_counter_ = nullptr;
		}
	}

	template<typename Type>
	Mat_<Type>::Mat_(int rows, int cols, int type) {
		create(rows, cols, type);
	}

	template<typename Type>
	Mat_<Type>::Mat_(Size size, int type)
		: Mat_(size.width_, size.height_, type) {}


 	template<typename Type>
 	Mat_<Type>::Mat_(int rows, int cols, int type, const Scalar& s) {
 		create(rows, cols, type);
 		*this = s;
 	}

 	template<typename Type>
 	Mat_<Type>::Mat_(Size size, int type, const Scalar& s) {
 		create(size.width_, size.height_, type);
 		*this = s;
 	}


	template<typename Type>
	Mat_<Type>::Mat_(const Mat_& m) {
		*this = m;
	}


	template<typename Type>
	Mat_<Type>::Mat_(int rows, int cols, int type, Type value)
		: Mat_(rows, cols, type) {
		set_value(value);
	}
	
	template<typename Type>
	Mat_<Type>& Mat_<Type>::operator=(const Mat_& m) {
		if (&m != this) {
			// 如果将m赋值给一当前的Mat（即m指向当前Mat），而且新的Mat的reference_counter_不为0，即
			// m还指向另一个Mat时（即m指向2个及以上的Mat），只能将m的reference counter加1；但如果m只指向当前
			// Mat，则直接将m的数据区直接销毁，并将m的各个attribute拷贝到当前Mat
			if (m.reference_counter_) {
				add_reference_counter(m.reference_counter_, 1);
			}
			// 如果新的Mat的reference_counter_为0，则先直接释放当前Mat
			//free();
			this->flags_ = m.flags_;
			this->rows_ = m.rows_;
			this->cols_ = m.cols_;
			this->layer_ = m.layer_;
			this->channels_ = m.channels_;
			this->channel_size_ = m.channel_size_;
			this->size_ = m.size_;
			this->data_ptr_ = m.data_ptr_;
			this->reference_counter_ = m.reference_counter_;
			this->data_begin_ = m.data_begin_;
			this->data_end_ = m.data_end_;
		} 
		return *this;
	}


	template<typename Type>
	Mat_<Type>::Mat_(const Mat_& m, const std::vector<int>& row, const std::vector<int>& col) {
		if (m.dims_ > 2) {
			std::vector<int> dimensions(m.dims_);
			dimensions[0] = row.begin();
			dimensions[1] = col.begin();
			for (int i = 2; i < m.dims_; i++) {
				//dimensions[i]=
			}
			*this = m(dimensions.data());
			return;
		}

		*this = m;

		if (!row.empty()) {
			this->rows_ = row.size();
			this->data_ptr_ += this->layer_*row.begin();
			this->flags_ |= (1 << 12);
		}

		if (!col.empty()) {
			this->cols_ = row.size();
			this->data_ptr_ += channel_size()*col.begin();
			this->flags_ |= (1 << 12);
		}
	}


// 	template<typename Type>
// 	Mat_<Type>::Mat_(int dims, const uint32_t* arr, int type) {
// 		assert(dims >= 0 && arr);
// 		// 如果试图创建1维Mat，会将其自动转换为2维，此时其第二维长度为1
// 		if ((dims == this->dims_ || (dims == 1 && dims_ <= 2)) && this->type_ = type()) {
// 
// 		}
// 	}
 	
 	template<typename Type>
  	Mat_<Type>& Mat_<Type>::operator=(const Scalar& s) {
 		if (this->empty())
 			return *this;
 		// FIXME: 此处未考虑多维平面
		size_t total_size = this->size_;		
		for (size_t i = 0; i < total_size; i++) {
			reinterpret_cast<uint8_t*>(data_ptr_)[i] = s.arr[i % 3];
		}
		return *this;
  	}

	template<typename Type>
	Mat_<Type>& Mat_<Type>::operator=(const std::initializer_list<int8_t>& value) {
		size_t index = 0;
		for (auto& iter : value) {
			this->data_ptr_[index] = iter;
			index++;
		}
		return *this;
	}

	template<typename Type>
	Mat_<Type>::~Mat_() {
		//free();
	}

	template<typename Type>
	int Mat_<Type>::channels() const {
		return this->channels_;
	}

	template<typename Type>
	int Mat_<Type>::type() const {
		return ((this->channels()) << 3 | this->depth());
	}

	template<typename Type>
	int Mat_<Type>::depth() const {
		return (this->flags_) & 0x07;
	}

	template<typename Type>
	size_t Mat_<Type>::element_size() const {
		return this->element_size_;
	}

	template<typename Type>
	size_t Mat_<Type>::channel_size() const {
		return this->channel_size_;
	}

	template<typename Type>
	bool Mat_<Type>::is_continuous() const {
		return ((this->flags_) >> 11) != 0;
	}

	template<typename Type>
	bool Mat_<Type>::is_submat() const {
		return ((this->flags_) >> 12) != 0;
	}

	template<typename Type>
	int Mat_<Type>::rows() const {
		return this->rows_;
	}

	template<typename Type>
	int Mat_<Type>::cols() const {
		return this->cols_;
	}

	template<typename Type>
	Mat_<Type> Mat_<Type>::col(int index) const {
		
	}


	template<typename Type>
	size_t Mat_<Type>::total() const {
		if (this->dims_ <= 2) {
			return this->size_;
		}
		size_t ret = 1;
		for (int i = 0; i < this->dims_; i++) {
			ret *= layer_[i];
		}
		return ret;
	}

	template<typename Type>
	size_t Mat_<Type>::total(int dim_begin, int dim_end) const {
		assert(dim_begin >= 0 && dim_beign <= dim_end);
		size_t ret = 1;
		int dim_end = (dim_end <= this->dims_ ? dim_end : this->dims_);
		for (int i = dim_end; i < dim_end; i++) {
			ret *= this->layer_[i];
		}
		return ret;
	}

	template<typename Type>
	size_t Mat_<Type>::nums() const {
		return total();
	}
	
	template<typename Type>
	Type& Mat_<Type>::at(int row, int col) {
		return *reinterpret_cast<Type*>(data_ptr_ + layer_[0] * row + layer_[1] * col);
	}

 	template<typename Type>
 	template<typename Type2>
 	Type2& Mat_<Type>::at(int row, int col) {
 		assert((row <= this->rows_) && (col <= this->cols_));
 		return *reinterpret_cast<Type2*>(data_ptr_ + layer_[0] * row + layer_[1] * col);
	}
 
  	template<typename Type>
  	template<typename Type2>
  	const Type2& Mat_<Type>::at(int row, int col) const {
  		assert((row <= this->rows_) && (col <= this->cols_));
  		return *reinterpret_cast<Type2*>(data_ptr_ + layer_[0] * row + layer_[1] * col);
	}


	template<typename Type>
	template<typename Type2>
	const Type2& Mat_<Type>::at(int row) const {
		assert(row <= this->rows_);
		return *reinterpret_cast<Type2>(data_ptr_ + layer_[0] * row);
	}

	template<typename Type>
	template<typename Type2>
	Type2& Mat_<Type>::at(int row) {
		assert(row <= this->rows_);
		return *reinterpret_cast<Type2>(data_ptr_ + layer_[0] * row);
	}

	template<typename Type>
	template<typename Type2>
	const Type2& Mat_<Type>::at(Point pt) const {
		return at<Type2>(pt.x_, pt.y_);
	}

	template<typename Type>
	template<typename Type2>
	Type2& Mat_<Type>::at(Point pt) {
		return at<Type>(pt.x_, pt.y_);
	}

// 	template<typename Type>
//  	template<int n>
//  	const Type& Mat_<Type>::at(Vec_<int, n>& vec) const {
//  		assert(sizeof(Type) == element_size());
// 		std::vector<uint32_t>& rgb;
// 		for (int i = 0; i < this->channels_; i++) {
// 			rgb.push_back()
// 		}
//  		//return *(Type*)(+vec.)
//  	}

	template<typename Type>
	bool Mat_<Type>::empty() const {
		return total() == 0;
	}

	template<typename Type>
	Mat_<Type> Mat_<Type>::zeros(int rows, int cols, int type) {
		return Mat_<Type>(rows, cols, type, 0);
	}

	template<typename Type>
	Mat_<Type> Mat_<Type>::zeros(Size size,int type) {
		return zeros(size.width_, size.height_, type);
	}

	template<typename Type>
	Mat_<Type> Mat_<Type>::ones(int rows, int cols, int type) {
		return Mat_<Type>(rows, cols, type, 1);
	}

	template<typename Type>
	Mat_<Type> Mat_<Type>::ones(Size size, int type) {
		return ones(size.width_, size.height_, type);
	}
	
	template<typename Type>
	void Mat_<Type>::set_value(Type value) {
		// 当有多行时
		const int total_size = nums();
			for (int i = 0; i < total_size; i++) {
				reinterpret_cast<Type*>(data_ptr_)[i] = value;
		}
	}


	template<typename Type>
	const Type* Mat_<Type>::operator[](size_t n) const {
		assert(n < this->rows_);
		return reinterpret_cast<Type*>(data_ptr_ + n * layer_[0]);
	}

	template<typename Type>
	Type* Mat_<Type>::operator[](size_t n) {
		assert(n < this->rows_);
		return reinterpret_cast<Type*>(data_ptr_ + n * layer_[0]);
	}


// 	template<typename Type>
// 	const Mat_<Type>& Mat_<Type>::clone() const {
// 		return Mat_<Type>(this->rows(), this->cols(), this->type());
// 	}

	template<typename Type>
	Mat_<Type>& operator-(Mat_<Type>& src, Mat_<Type>& dest) {
		//Mat_<Type> tmp(src.clone());
		int channels = src.channels();
		for (int i = 0; i < src.rows(); i++) {
			//const uint8_t* src_row = src.ptr(i);
			for (int j = 0; j < src.cols(); j++) {
				for (int c = 0; c < channels; c++) {
					dest.at<Vec3uc>(i, j)[c] -= src.at<Vec3uc>(i, j)[c];
				}
			}
		}
		return dest;
	}

	template<typename Type>
	const Size Mat_<Type>::size() const {
		return Size(this->rows(), this->cols());
	}

	template<typename Type>
	Type* Mat_<Type>::ptr(int row) {
		assert(row < this->rows_);
		return reinterpret_cast<Type*>(data_ptr_ + row * layer_[0]);
	}

	template<typename Type>
	const Type* Mat_<Type>::ptr(int row) const {
		assert(row < this->rows_);
		return reinterpret_cast<Type*>(data_ptr_ + row * layer_[0]);
	}

	template<typename Type>
	template<typename Type2>
	Type2* Mat_<Type>::ptr(int row) {
		assert(row < this->rows_);
		return reinterpret_cast<Type2*>(data_ptr_ + row * layer_[0]);
	}

	template<typename Type>
	template<typename Type2>
	const Type2* Mat_<Type>::ptr(int row) const {
		assert(row < this->rows_);
		return reinterpret_cast<Type2*>(data_ptr_ + row * layer_[0]);
	}

	


	//////////////////// mat_const_iterator //////////////
 	template<typename Type>
	mat_const_iterator<Type>::mat_const_iterator()
		: start_(0), end_(0), cur_(0), mat_(0), element_size_(0) {}
	
	template<typename Type>
	mat_const_iterator<Type>::mat_const_iterator(const Mat_<Type>& m)
		: mat_(m), element_size_(m.element_size()) {
		if (m && m.is_continuous()) {
			this->begin_ = m.data_ptr_;
			this->end_ = this->begin_ + m.size();
		}
	}

	template<typename Type>
	mat_const_iterator<Type>::mat_const_iterator(const mat_const_iterator<Type>& other)
		: start_(iter.begin_), end_(iter.end_), cur_(iter.cur_), mat_(iter.mat_), element_size_(iter.element_size_) {}

	template<typename Type>
	mat_const_iterator<Type>& mat_const_iterator<Type>::operator=(const mat_const_iterator<Type>& other) {
		this->begin_ = other.begin_;
		this->end_ = other.end_;
		this->mat_ = other.mat_;
		this->cur_ = other.cur_;
		this->element_size_ = other.element_size_;
		return *this;
	}

	template<typename Type>
	const Type& mat_const_iterator<Type>::operator*() const {
		return this->cur_;
	}

	template<typename Type>
	const Type& mat_const_iterator<Type>::operator[](difference_type diff) const {
		return this->cur_ + diff;
	}

// 	template<typename Type>
// 	const iterator mat_const_iterator<Type>::begin() const {
// 		return this->begin_;
// 	}
// 
// 	template<typename Type>
// 	const iterator mat_const_iterator<Type>::end() const {
// 		return this->end_;
// 	}

}
