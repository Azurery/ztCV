#include "Traits.h"
#include <utility>
#include <cassert>

namespace ztCV {

	template<typename Type> Point_<Type>::Point_()
		: x_(0), y_(0) {}

	template<typename Type> Point_<Type>::Point_(Type x, Type y)
		: x_(x), y_(y) {}


	template<typename Type> Point_<Type>::Point_(const Point_& pt)
		: x_(std::move(pt.x_)), y_(std::move(pt.y_)) {}

	template<typename Type> Point_<Type>::Point_(Point_&& pt) noexcept
		: Point_(std::move(pt)) {}

	template<typename Type> Point_<Type>& Point_<Type>::operator=(const Point_& pt) {
		x_ = pt.x_;
		y_ = pt.y_;
		return *this;
	}

	template<typename Type> Point_<Type>& Point_<Type>::operator=(const Point_&& pt) {
		x_ = std::move(pt.x_);
		y_ = std::move(pt.y_);
		return *this;
	}

	template<typename Type> Type Point_<Type>::inner_product(const Point_& pt) {
		return static_cast<Type>(this->x_)*(pt.x_) + static_cast<Type>(this->y_)*(pt.y_);
	}

	template<typename Type> double Point_<Type>::outer_product(const Point_& pt) {
		return static_cast<double>(this->x_)*(pt.y_) + static_cast<double>(this->y_)*(pt.x_);
	}

	template<typename Type> static 
	Point_<Type>& operator+=(Point_<Type>& pt1, const Point_<Type>& pt2) {
		pt1.x_ += pt2.x_;
		pt1.y_ += pt2.y_;
		return pt1;
	}

	template<typename Type> static
		Point_<Type>& operator-=(Point_<Type>& pt1, const Point_<Type>& pt2) {
		pt1.x -= pt2.x_;
		pt1.y_ -= pt2.y_;
		return pt1;
	}


	//////////////// Vec_ /////////////////////////
	template<typename Type, int n> Vec_<Type, n>::Vec_() {
// 		for (int i = 0; i < n; i++) {
// 			arr[i] = 0;
// 		}
		memset(arr, 0, n);
	}

	template<typename Type, int n> Vec_<Type, n>::Vec_(Type vec0) {
		assert(n >= 1);
		arr[0] = vec0;
// 		for (int i = 1; i < n; i++) {
// 			arr[i] = 0;
// 		}
		memset(arr + 1, 0, n - 1);
	}

	template<typename Type, int n> Vec_<Type, n>::Vec_(Type vec0, Type vec1) {
		static_assert(n >= 2);
		arr[0] = vec0, arr[1] = vec1;
// 		for (int i = 2; i < n; i++) {
// 			arr[i] = 0;
// 		}
		memset(arr + 2, 0, n - 2);
	}

	template<typename Type, int n> Vec_<Type, n>::Vec_(Type vec0, Type vec1, Type vec2) {
		assert(n >= 3);
		arr[0] = vec0, arr[1] = vec1, arr[2] = vec2;
		for (int i = 3; i < n; i++) {
			arr[i] = 0;
		}
	}

	template<typename Type, int n> Vec_<Type, n>::Vec_(Type vec0, Type vec1, Type vec2, Type vec3) { 
		static_assert(n >= 4);
		arr[0] = vec0, arr[1] = vec1, arr[2] = vec2, arr[3] = vec3;
		for (int i = 4; i < n; i++) {
			arr[i] = 0;
		}
	}

	template<typename Type, int n> Vec_<Type, n>::Vec_(Vec_<Type, n>& vec) 
		: Vec_<Type, n>(vec.arr[0], vec.arr[1], vec.arr[2]) {}

	template<typename Type, int n> Vec_<Type, n>::Vec_(std::initializer_list<Type>& il) {
		int i = 0;
		for (auto iter : il) {
			arr[i++] = iter;
		}
	}

	template<typename Type, int n>
	const Type& Vec_<Type, n>::operator[](int index) const {
		assert(index < n);
		return arr[index];
	}
	
	template<typename Type, int n>
	Type& Vec_<Type, n>::operator[](int index) {
		assert(index < n);
		return arr[index];
	}

	template<typename Type, int n>
	const Type& Vec_<Type, n>::operator()(int index) const {
		static_assert(index < n);
		return arr[i];
	}

	template<typename Type, int n>
	Type& Vec_<Type, n>::operator()(int index) {
		static_assert(index < n);
		return arr[i];
	}

	template<typename Type>
	Size_<Type>::Size_()
		: width_(0), height_(0) {}

	template<typename Type>
	Size_<Type>::Size_(Type width, Type height)
		: width_(width), height_(height) {}

	template<typename Type>
	Size_<Type>::Size_(const Size_& s)
		: Size_(s.width_, s.height_) {}

	template<typename Type>
	Size_<Type>::Size_(Size_&& s) noexcept
		: Size_(std::move(s.width_), std::move(s.height_)) {}

	template<typename Type>
	Size_<Type>::Size_(const Point_<Type>& pt)
		: Size_(pt.x_, pt.y_) {}

	template<typename Type>
	Size_<Type>& Size_<Type>::operator=(const Size_& s) {
		width_ = s.width_;
		height_ = s.height_;
		return *this;
	}

	template<typename Type>
	Size_<Type>& Size_<Type>::operator=(Size_&& s) noexcept {
		width_ = std::move(s.width_);
		height_ = std::move(s.height_);
		return *this;
	}

	template<typename Type>
	bool Size_<Type>::operator==(const Size_& s) {
		return (this->height_ == s.height_  && this->width_ == s.width_);
	}

	template<typename Type>
	Type Size_<Type>::area() const {
		return width_ * height_;
	}

	template<typename Type>
	bool Size_<Type>::empty() const {
		return (width_ <= 0 || height_ <= 0);
	}
	
	template<typename Type> template<typename Type2>
	Size_<Type>::operator Size_<Type2>() const {
		return Size_<Type2>(width_, height_);
	}

	template<typename Type>
	Rect_<Type>::Rect_()
		: x_(0), y_(0), width_(0), height_(0) {}

	template<typename Type>
	Rect_<Type>::Rect_(Type x, Type y, Type width, Type height)
		: x_(x), y_(y), width_(width), height_(height) {}

	template<typename Type>
	Rect_<Type>::Rect_(const Rect_& r)
		: Rect_(r.width_, r.height_) {}

	template<typename Type>
	Rect_<Type>::Rect_(Rect_&& r) noexcept
		: Rect_(std::move(r.width_), std::move(r.height_)) {}

	template<typename Type>
	Rect_<Type>::Rect_(const Point_<Type>& pt1, const Point_<Type>& pt2) {
		this->x_ = pt1.x_;
		this->y_ = pt1.y_;
		this->width_ = pt2.x_ - pt1.x_;
		this->height_ = pt2.y_ - pt2.y_;
	}

	template<typename Type>
	Rect_<Type>& Rect_<Type>::operator=(const Rect_& r) {
		x_ = r.x_;
		y = r.y_;
		return *this;
	}

	template<typename Type>
	Rect_<Type>& Rect_<Type>::operator=(Rect_&& r) noexcept {
		x_ = std::move(r.x_);
		y_ = std::move(r.y_);
		return *this;
	}

	template<typename Type> template<typename Type2>
	Rect_<Type>::operator Rect_<Type2>() const {
		return Rect_<Type2>(this->x_, this->y_);
	}

	template<typename Type>
	Scalar_<Type>::Scalar_() {
		for (int i = 0; i < 4; i++) {
			this->arr[i] = 0;
		}
	}

	template<typename Type>
	Scalar_<Type>::Scalar_(Type t0, Type t1, Type t2) {
		this->arr[0] = t0;
		this->arr[1] = t1;
		this->arr[2] = t2;
	}

	template<typename Type>
	Scalar_<Type>::Scalar_(const Scalar_& s) : Vec_<Type, 4>(s) {}

	template<typename Type>
	Scalar_<Type>::Scalar_(Scalar_&& s) noexcept : Vec_<Type, 4>(std::move(s)) {}

	template<typename Type>
	Scalar_<Type>& Scalar_<Type>::operator=(const Scalar_& s) {
		for (int i = 0; i < 4; i++) {
			this->arr[i] = s.arr[i];
		}
		return *this;
	}

	template<typename Type>
	Scalar_<Type>& Scalar_<Type>::operator=(Scalar_&& s) noexcept {
		for (int i = 0; i < 4; i++) {
			this->arr[i] = std::move(s.arr[i]);
		}
		return *this;
	}

	template<typename Type> template<typename Type2>
	Scalar_<Type>::operator Scalar_<Type2>() const {
		return Scalar_<Type2>(this);
	}

	template<typename Type>
	Complex_<Type>::Complex_()
		: real_(0), imaginary_(0) {}

	template<typename Type>
	Complex_<Type>::Complex_(Type real, Type imaginary)
		: real_(real), imaginary_(imaginary) {}

	template<typename Type>
	Complex_<Type> Complex_<Type>::conjugate() const {
		return Complex_<Type>(this->real_, -(this->imaginary_));
	}


}
