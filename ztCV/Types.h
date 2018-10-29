#ifndef _ZT_CV_TYPES_H_
#define _ZT_CV_TYPES_H_
#include <initializer_list>
#include <cstdlib>
#include <cstdint>
#include "Traits.h"
namespace ztCV {
	//////////////////////Point_/////////////////////////
	/*
	@brief Point_���ڱ�ʾ��ά�ĵ㣬ӵ��x��y����

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
	@brief  Vec_��Ҫ���ڶ�ͨ��ͼ��

	����Vec����࣬����������ʾ����������������ֵ���飬���ҿ��Զ��������ѧ���㣬
	��ȡԪ�صȡ�

	��Ҫ��������������
	@param Type �� Ԫ������
	@param n    �� Ԫ������

	ͨ������ʹ�÷�ʽΪ Vec<int,3> ������ӵ��3��int���͵�Ԫ�ص����飬�����Ը��Ӽ�
	��д��Vec3i��������ȫ�ȼۡ�

	֧�ֵĲ�����
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
	* @brief Size_����ָ��һ��image����rectangle��size
	*
	* @description ���class����������member��width��height
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
	* @brief Rect_���class���ڱ�ʾ���Σ�rectangle��
	*
	* @description :
	*		1. width��height
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
	* @brief Scalar_���classӵ��4��Ԫ�أ��̳���Vec
	*
	* @description ���class�̳��� Vec<Type, 4> ��Scalar_��Scalar����
	*				��������һ�����͵�vector��ӵ��4��Ԫ�ء�Ĭ��ֵΪ0������������ֵ��
	*				������ͨ����ֵ����һ�����壬����������ͬͨ�����ľ���Ԫ�ء�
	* @example Mat M(7,7,CV_32FC2,Scalar(1,3));
	*			����һ��2ͨ������ÿ��ͨ����ֵ��Ϊ��1,3�������Ϊ32��7��7�е�ͼ�����
	*				CV_32F��ʾÿ��Ԫ�ص�ֵ������Ϊ32λ��������C2��ʾͨ����Ϊ2��
	*			Scalar��1,3����ʾ�Ծ���ÿ��Ԫ�ض���ֵΪ��1,3������һ��ͨ���е�ֵ����1��
	*				�ڶ���ͨ���е�ֵ����3.	
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
	* @brief Complex_���class���ڱ�ʾ����
	*
	* @description ���class��std::complex���Ʋ��Ҽ��ݣ��������Ը��ӷ���ط���һ��������
	*				��ʵ����real�����鲿��imaginary��
	*
	*/
	template<typename Type>
	class Complex_ {
	public:
		Complex_();
		Complex_(Type real, Type imaginary = 0);

		template<typename Type2> operator Complex_<Type2>() const;
		
		// ���ڷ��ع����
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