#ifndef _ZT_CV_TYPES_H_
#define _ZT_CV_TYPES_H_
#include <initializer_list>
#include <cstdlib>
#include <cstdint>
#include "Traits.h"
namespace ztCV {
	//////////////////////Point_/////////////////////////
	/*
	@brief Point_用于表示二维的点，拥有x和y坐标

	@operator implementation
	point1 = point2 + point3
	point1 = point2 - point3
	point1 = num * point2
	point1 = point2 * num
	point1 = point2 / num
	point1 += point2
	point1 -= point2
	point1 *= point2
	point1 /= point2
	point1 == point2
	point1 != point2

	@type aliases
	typedef Point_<int> Point2i;
	typedef Point2i Point;
	typedef Point_<float> Point2f;
	typedef Point_<double> Point2d;
	*/
	template<typename Type>
	class Point_ {
	public:
		using value_type = Type;

		Point_();
		Point_(Type x, Type y);
		Point_(const Point_& pt);
		Point_(Point_&& pt) noexcept;

		Point_& operator=(const Point_& pt);
		Point_& operator=(const Point_&& pt);

		Type inner_product(const Point_& pt);
		double outer_product(const Point_& pt);
		Type x_;
		Type y_;
	};

	using Point2i = Point_<int>;
	using Point2d = Point_<double>;
	using Point2f = Point_<float>;
	using Point = Point2i;

	//////////////////////Vec_////////////////
	/*
	@brief  Vec_主要用于多通道图像

	利用Vec这个类，可以用来表示具有有限数量的数值数组，并且可以对其进行数学运算，
	获取元素等。

	需要传入两个参数：
	@param Type ： 元素类型
	@param n    ： 元素数量

	通常，其使用方式为 Vec<int,3> ，代表拥有3个int类型的元素的数组，但可以更加简化
	地写作Vec3i，两者完全等价。

	支持的操作：
	-   v1 = v2 + v3
	-   v1 = v2 - v3
	-   v1 = v2 \* scale
	-   v1 = scale \* v2
	-   v1 = -v2
	-   v1 += v2 and other augmenting operations
	-   v1 == v2, v1 != v2
	*/
	template<typename Type, int n>
	class Vec_ {
	public:
		using value_type = Type;
		Vec_();
		Vec_(Type vec0);
		Vec_(Type vec0, Type vec1);
		Vec_(Type vec0, Type vec1, Type vec2);
		Vec_(Type vec0, Type vec1, Type vec2, Type vec3);

		Vec_(Vec_<Type, n>& vec);
		Vec_(std::initializer_list<Type>&);

		const Type& operator[](int index) const;
		Type& operator[](int index);
		const Type& operator()(int index) const;
		Type& operator()(int index);
		Type arr[n];
	};

	/// unsigned char ====> 8 bits
	using Vec2uc = Vec_<uint8_t, 2>;
	using Vec3uc = Vec_<uint8_t, 3>;
	using Vec4uc = Vec_<uint8_t, 4>;

	/// signed char  =====> 8 bits
	using Vec2sc = Vec_<int8_t, 2>;
	using Vec3sc = Vec_<int8_t, 3>;
	using Vec4sc = Vec_<int8_t, 4>;

	/// unsigned short int =====> 16 bits
	using Vec2us = Vec_<uint16_t, 2>;
	using Vec3us = Vec_<uint16_t, 3>;
	using Vec4us = Vec_<uint16_t, 4>;

	/// signed short int =====>16 bits
	using Vec2ss = Vec_<int16_t, 2>;
	using Vec3ss = Vec_<int16_t, 3>;
	using Vec4ss = Vec_<int16_t, 4>;

	/// int ======> 32 bits
	using Vec2i = Vec_<int32_t, 2>;
	using Vec3i = Vec_<int32_t, 3>;
	using Vec4i = Vec_<int32_t, 4>;
	using Vec6i = Vec_<int32_t, 6>;
	using Vec8i = Vec_<int32_t, 8>;

	/// float =====> 32 bits
	using Vec2f = Vec_<float, 2>;
	using Vec3f = Vec_<float, 3>;
	using Vec4f = Vec_<float, 4>;
	using Vec6f = Vec_<float, 6>;

	/// double ======> 64 bits
	using Vec2d = Vec_<double, 2>;
	using Vec3d = Vec_<double, 3>;
	using Vec4d = Vec_<double, 4>;
	using Vec6d = Vec_<double, 6>;
	
	///////////////// Size_ /////////////////////
	/**
	* @brief Size_用于指定一个image或者rectangle的size
	*
	* @description 这个class包含了两个member：width和height
	*
	*/
	template<typename Type>
	class Size_ {
	public:
		using value_type = Type;

		Size_();
		Size_(Type width, Type height);
		Size_(const Size_& s);
		Size_(Size_&& s) noexcept;
		Size_(const Point_<Type>& pt);

		Size_& operator=(const Size_& s);
		Size_& operator=(Size_&& s) noexcept;
		bool operator==(const Size_& s);

		Type area() const;
		bool empty() const;

		template<typename Type2> operator Size_<Type2>() const;

		Type width_;
		Type height_;
	};

	using Size2i = Size_<int>;
	using Size2f = Size_<float>;
	using Size2d = Size_<double>;
	using Size = Size2i;

	template<typename Type>
	class DataType<Size_<Type>> {
	public:
		using value_type = Size_<Type>;
		using channel_type = Type;
		enum {
			channels = 2,
			depth = DataType<channel_type>::depth,
		    type = MAKE_TYPE(depth, channels)
		};
	};

	////////////////////// Rect_ ///////////////////
	/**
	* @brief Rect_这个class用于表示矩形（rectangle）
	*
	* @description :
	*		1. width和height
	*              
	*
	*/
	template<typename Type>
	class Rect_ {
	public:
		using value_type = Type;
		Rect_();
		Rect_(Type x, Type y, Type width, Type height);

		Rect_(const Rect_& r);
		Rect_(Rect_&& r) noexcept;
		Rect_(const Point_<Type>& pt1, const Point_<Type>& pt2);

		Rect_& operator=(const Rect_& r);
		Rect_& operator=(Rect_&& r) noexcept;

		template<typename Type2> operator Rect_<Type2>() const;

		Type x_;
		Type y_;
		Type width_;
		Type height_;
	};


	template<typename Type>
	class DataType<Rect_<Type>> {
	public:
		using value_type = Rect_<Type>;
		using channel_type = Type;
		enum {
			channels = 4,
			depth = DataType<channel_type>::depth,
			type = MAKE_TYPE(depth, channels)
		};
	};

	//////////////// Scalar_ /////////////////
	/**
	* @brief Scalar_这个class拥有4个元素，继承自Vec
	*
	* @description 这个class继承自 Vec<Type, 4> ，Scalar_和Scalar可以
	*				被视作是一个典型的vector，拥有4个元素。默认值为0（如果不传入初值）
	*				将各个通道的值构成一个整体，赋给具有相同通道数的矩阵元素。
	* @example Mat M(7,7,CV_32FC2,Scalar(1,3));
	*			创建一个2通道，且每个通道的值都为（1,3），深度为32，7行7列的图像矩阵。
	*				CV_32F表示每个元素的值的类型为32位浮点数，C2表示通道数为2，
	*			Scalar（1,3）表示对矩阵每个元素都赋值为（1,3），第一个通道中的值都是1，
	*				第二个通道中的值都是3.	
	*/
	template<typename Type>
	class Scalar_ : public Vec_<Type, 4> {
	public:
		Scalar_();
		Scalar_(Type t0, Type t1, Type t2);
		Scalar_(const Scalar_& s);
		Scalar_(Scalar_&& s) noexcept;

		Scalar_& operator=(const Scalar_& s);
		Scalar_& operator=(Scalar_&& s) noexcept;

		template<typename Type2> operator Scalar_<Type2>() const;

		static Scalar_<Type> set_all(Type t);
	};

	using Scalar = Scalar_<uint32_t>;

	template<typename Type>
	class DataType<Scalar_<Type>> {
	public:
		using value_type = Scalar_<Type>;
		using channel_type = Type;
		enum {
			channels = 4,
			depth = DataType<channel_type>::depth,
			type = MAKE_TYPE(depth, channels)
		};
	};


	//////////////////// Complex_ //////////////////////////
	/**
	* @brief Complex_这个class用于表示复数
	*
	* @description 这个class与std::complex相似并且兼容，但它可以更加方便地访问一个复数的
	*				的实部（real）和虚部（imaginary）
	*
	*/
	template<typename Type>
	class Complex_ {
	public:
		Complex_();
		Complex_(Type real, Type imaginary = 0);

		template<typename Type2> operator Complex_<Type2>() const;
		
		// 用于返回共轭复数
		// @example a+bi =====> a-bi
		Complex_ conjugate() const;

		Type real_;
		Type imaginary_;
	};

	using Complexf = Complex_<float>;
	using Complexd = Complex_<double>;

	template<typename Type>
	class DataType<Complex_<Type>> {
	public:
		using value_type = Complex_<Type>;
		using channel_type = Type;
		enum {
			channels = 2,
			depth = DataType<channel_type>::depth,
			type = MAKE_TYPE(depth, channels)
		};
	};

}
#endif
//#include "Types.hpp"