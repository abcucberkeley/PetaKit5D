# include <mex.h>

# include <image.hpp>

// Matlab mxArray size are row-wise. This trait provides the following
// conversion:
// [nrows, ncols] -> [width, height]
// [nrows, ncols, nslices] -> [width, height, depth]

template <int n>
struct sizeWrapper
{
  template <typename S1, typename S2>
  static bool eq(const S1* s1, const S2* s2);
  
  template <typename S1, typename S2>
  static void convert(const S1* src, S2* dst);
};
 
template <>
template <typename S1, typename S2>
void sizeWrapper<2>::convert(const S1* src, S2* dst)
{
  dst[0] = src[1];
  dst[1] = src[0];
}

template <>
template <typename S1, typename S2>
bool sizeWrapper<2>::eq(const S1* s1, const S2* s2)
{
  return s1[0] == s2[1] && s1[1] == s2[0];
}

template <>
template <typename S1, typename S2>
void sizeWrapper<3>::convert(const S1* src, S2* dst)
{
  dst[0] = src[1];
  dst[1] = src[0];
  dst[2] = src[2];
}

template <>
template <typename S1, typename S2>
bool sizeWrapper<3>::eq(const S1* s1, const S2* s2)
{
  return s1[0] == s2[1] && s1[1] == s2[0] && s1[2] == s2[2];
}

// This trait makes the correspondance between C built-in types
// (e.g. int, double, etc.) amd Matlab type (e.g. mxDOUBLE_CLASS).

template <typename T>
struct mxType
{
  static const mxClassID mxID;
};

template <> const mxClassID mxType<char>::mxID = mxINT8_CLASS;
template <> const mxClassID mxType<unsigned char>::mxID = mxUINT8_CLASS;
template <> const mxClassID mxType<short>::mxID = mxINT16_CLASS;
template <> const mxClassID mxType<unsigned short>::mxID = mxUINT16_CLASS;
template <> const mxClassID mxType<int>::mxID = mxINT32_CLASS;
template <> const mxClassID mxType<unsigned int>::mxID = mxUINT32_CLASS;
template <> const mxClassID mxType<long>::mxID = mxINT64_CLASS;
template <> const mxClassID mxType<unsigned long>::mxID = mxUINT64_CLASS;
template <> const mxClassID mxType<float>::mxID = mxSINGLE_CLASS;
template <> const mxClassID mxType<double>::mxID = mxDOUBLE_CLASS;

// This function copy the content of an image object into a mxArray
// from Matlab.

template <int n, typename T>
void image2mxArray(const image<n, T> & src, mxArray*& dst)
{
  const int* src_size = src.size();
  mwSize src_size2[n];
  sizeWrapper<n>::convert(src_size, src_size2);

  if (dst)
    {
      if (mxGetNumberOfDimensions(dst) != n ||
	  !sizeWrapper<n>::eq(src_size, mxGetDimensions(dst)))
	{
	  mxDestroyArray(dst);
	  dst = mxCreateNumericArray(n, src_size2, mxType<T>::mxID, mxREAL);
	}
    }
  else
    dst = mxCreateNumericArray(n, src_size2, mxType<T>::mxID, mxREAL);
  
  src.raw_data((T*) mxGetData(dst));
}
