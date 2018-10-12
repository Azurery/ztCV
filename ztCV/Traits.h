#ifndef _ZT_CV_TRAITS_H_
#define _ZT_CV_TRAITS_H_
#include <cstdint>
#include <cstdlib>


#define MAT_FLAGS_MASK (1<<12)
#define MAKE_TYPE(depth,channels) ((depth)&(00000111)|((channels)<<3))

#define CV_8UCC1 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT8),1))
#define CV_8UCC2 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT8),2))
#define CV_8UCC3 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT8),3))
#define CV_8UCC4 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT8),4))

#define CV_8SCC1 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT8),1))
#define CV_8SCC2 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT8),2))
#define CV_8SCC3 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT8),3))
#define CV_8SCC4 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT8),4))

#define CV_16USC1 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT16),1))
#define CV_16USC2 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT16),2))
#define CV_16USC3 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT16),3))
#define CV_16USC4 (MAKE_TYPE(static_cast<int>(cv_depth::CV_UINT16),4))

#define CV_16SSC1 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT16),1))
#define CV_16SSC2 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT16),2))
#define CV_16SSC3 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT16),3))
#define CV_16SSC4 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT16),4))

#define CV_32SIC1 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT32),1))
#define CV_32SIC2 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT32),2))
#define CV_32SIC3 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT32),3))
#define CV_32SIC4 (MAKE_TYPE(static_cast<int>(cv_depth::CV_INT32),4))

#define CV_32FC1 (MAKE_TYPE(static_cast<int>(cv_depth::CV_FLOAT),1))
#define CV_32FC2 (MAKE_TYPE(static_cast<int>(cv_depth::CV_FLOAT),2))
#define CV_32FC3 (MAKE_TYPE(static_cast<int>(cv_depth::CV_FLOAT),3))
#define CV_32FC4 (MAKE_TYPE(static_cast<int>(cv_depth::CV_FLOAT),4))

#define CV_64DC1 (MAKE_TYPE(static_cast<int>(cv_depth::CV_DOUBLE),1))
#define CV_64DC2 (MAKE_TYPE(static_cast<int>(cv_depth::CV_DOUBLE),2))
#define CV_64DC3 (MAKE_TYPE(static_cast<int>(cv_depth::CV_DOUBLE),3))
#define CV_64DC4 (MAKE_TYPE(static_cast<int>(cv_depth::CV_DOUBLE),4))

namespace ztCV {
	/**
	*  @brief  cv_depth这个enum表示图像深度。
	*
	*  @intepretion  图像深度：用于表示存储每个像素所用的位数，也用于
	*              度量图像色彩的分辨率。
	*                  图像深度确定彩色图像的每个像素可能有的颜色数，或者确定灰度图像的每个像素
	*              可能有的灰度级数，它决定了彩色图像中可能出现的最多的颜色数，或者灰度图像中的
	*                  最大灰度等级。
	* @example
	*
	*
	*/
	enum class cv_depth {
		CV_DEFAULT = -1,
		CV_UINT8 = 1,
		CV_INT8 = 1,
		CV_UINT16 = 2,
		CV_INT16 = 2,
		CV_INT32 = 4,
		CV_FLOAT = 4,
		CV_DOUBLE = 8,
	};

	/**
	* @brief DataType这个class其实就是一个C++基础类型（primitive）的包装类
	*
	* @description 定义DataType的目的是为了防止用户滥用数据，对于每一种C++基础数据
	*                  类型都有与之对应的DataType数据类型，而在DataType数据类型中
	*              也限定了几个基础类型（value_type chennel_type channel type depth），
	*                  即，一旦传入一个C++的基础数据类型，则在CV库中就会将与之对应的
	*              库类型与之绑定，从而防止了类型之间相互污染，
	*
	*  @example 如果想要创建一个浮点型的矩阵：
	*              Mat A(30,40,DataType<float>::type);
	*          一旦确定了矩阵的类型，该矩阵的各个参数便都确定了下来
	*              不存在用户去随意指定该矩形的各个参数
	*          其实，DataType是CV库与C++之间的一个桥梁
	*/
	template<typename Type>
	class DataType {
	public:
		using value_type = Type;
		using channel_type = Type;
		enum {
			depth = cv_depth::CV_DEFAULT,
			channels = -1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<bool> {
	public:
		using value_type = bool;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_UINT8,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<uint8_t> {
	public:
		using value_type = uint8_t;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_UINT8,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<int8_t> {
	public:
		using value_type = int8_t;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_UINT8,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<char> {
	public:
		using value_type = char;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_UINT8,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<uint16_t> {
	public:
		using value_type = uint16_t;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_UINT16,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<int16_t> {
	public:
		using value_type = int16_t;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_INT16,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<int32_t> {
		using value_type = int32_t;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_INT32,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<float> {
	public:
		using value_type = float;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_FLOAT,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};

	template<>
	class DataType<double> {
	public:
		using value_type = double;
		using channel_type = value_type;
		enum {
			depth = cv_depth::CV_DOUBLE,
			channels = 1,
			type = MAKE_TYPE(static_cast<int>(depth), channels)
		};
	};
}
#endif
