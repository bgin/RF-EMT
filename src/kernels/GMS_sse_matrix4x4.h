#ifndef __GMS_SSE_MATRIX_4X4_H__
#define __GMS_SSE_MATRIX_4X4_H__

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Generic templated Matrix 4x4  classes based upon SIMD types __m256d and __m256 for adequate specializations - classes declaration.

@file SIMDMatrix4x4.h
@author: Bernard Gingold
@version:  1.0  25/12/2015
@description: SIMDMatrix4x4.h
*/




namespace  gms 
{

namespace math
{

	/*
	     ***   Primary template  ***
	*/
	 
	template<class Args, typename T>  class SIMDMat4x4;
	
	
	

	/*
	      ***  Partial Float specialization  ***
	*/
	template<class Args>  class SIMDMat4x4<Args,float>
	{


	public:

		/*
		@brief   Default Ctor constructs , zero-filled matrix.
		*/
		SIMDMat4x4();

			
		

		/*
		@brief   Ctor constructs an object from the user passed 4x4 scalars wrapped 
		          in auxilary structure.
		*/
		SIMDMat4x4(_In_   const   Args &);

		/*
		@brief    Ctor constructs an object , from the user passed __m128 vector type.
		*/
		SIMDMat4x4(_In_ __m128 const&);

		/*
		@brief    Ctor constructs an object  from the user passed std::initializer_list
		*/
		SIMDMat4x4(_In_ std::initializer_list<__m128> const &);

		/*
		@brief    Ctor constructs an object, form the user passed scalar.
		*/
		SIMDMat4x4(_In_ const float);

		/*
		@brief    Ctor constructs an object form the user passed reference 
		          to static array of type float.
		*/
		SIMDMat4x4(_In_ float(&)[16]);

		/*
		@brief    Ctor constructs an object from the user passed static array
		          of type __m128 [4].
		*/
		SIMDMat4x4(_In_ __m128(&)[4]);

		/*
		@brief    Copy-Ctor.
		*/
		SIMDMat4x4(_In_ SIMDMat4x4 const&);

		/*
		@brief    Move-Ctor.
		*/
		SIMDMat4x4(_In_ SIMDMat4x4 &&);


		/*
		        ***    Class friend and member operators   ***
		*/

		/*
		@brief    operator= copy-assign.
		*/
		SIMDMat4x4 &         operator=(_In_ SIMDMat4x4 const &);

		/*
		@brief    operator=  move-assign.
		*/
		SIMDMat4x4 &         operator=(_In_ SIMDMat4x4 &&);

		/*
		@brief   operator+= , type arg const SIMDMat4x4 &
		*/
		SIMDMat4x4 &         operator+=(_In_ SIMDMat4x4 const&);

		/*
		@brief   operator+= , type arg const __m128 & .
		*/
		SIMDMat4x4 &         operator+=(_In_ __m128 const &);

		/*
		@brief   operator+= , type arg single scalar.
		*/
		SIMDMat4x4 &         operator+=(_In_ const float);

		/*
		@brief   operator-= , type arg const SIMDMat4x4 &.
		*/
		SIMDMat4x4 &         operator-=(_In_ SIMDMat4x4 const &);

		/*
		@brief   operator-= , type arg const __m128 &,
		*/
		SIMDMat4x4 &         operator-=(_In_ const __m128 &);

		/*
		@brief   operator-= , type arg is const float
		*/
		SIMDMat4x4 &         operator-=(_In_ const float);

		/*
		@brief   operator*= , type arg is const SIMDMat4x4 &. // ERROR WRONG IMPLEMENTATION.
		*/
		SIMDMat4x4 &         operator*=(_In_  SIMDMat4x4 &);

		/*
		@brief   operator*= , type arg is const __m128 &.
		*/
		SIMDMat4x4 &         operator*=(_In_ const __m128 &);

		/*
		@brief   operator*= , type arg is const float.
		*/
		SIMDMat4x4 &         operator*=(_In_ const float);

		/*
		@brief   operator/=, type arg is const float only.
		*/
		SIMDMat4x4 &         operator/=(_In_ const float);

		/*
		@brief   non-const operator[]().
		*/
		__m128               operator[](_In_ const int);

		/*
		@brief   const operator[]() const.
		*/
		const   __m128       operator[](_In_ const int) const;

		/*
		@brief    operator<<
		*/
		friend  std::ostream &      operator<<(_In_ std::ostream &, _In_ const SIMDMat4x4<Args,float> &);
		/*
		@brief   operator== , for arg type SIMDMat4x4 const &, based on SSE CMPPS instruction.
		*/
		SIMDMat4x4           operator==(_In_ SIMDMat4x4 const &);

		/*
		@brief   operator== ,for arg type __m128 const &, based on SSE CMPPS instruction.
		*/
		SIMDMat4x4           operator==(_In_ __m128 const &);

		/*
		@brief   operator== ,for arg type: const float, based on SSE CMPPS instruction.
		*/
		SIMDMat4x4           operator==(_In_ const float);

		/*
		@brief   operator!= ,for arg type: SIMDMat4x4 const &, based on SSE CMPPS instruction.
		*/
		SIMDMat4x4           operator!=(_In_ SIMDMat4x4 const &);

		/*
		@brief   operator!= , for arg type: __m128 const &, based on SSE CMPPS instruction.
		*/
		SIMDMat4x4           operator!=(_In_ __m128 const &);

		/*
		@brief   operator!= , for arg type: const float , based on SSE CMPPS instruction.
		*/
		SIMDMat4x4           operator!=(_In_ const float);

		/*
		         ***      Class member functions   ***
		*/

		/*
		@brief     inserts __m128 element at specified location.
		*/
		void      insert_at(_In_  int, _In_ __m128 &);

		/*
		@brief  Matrix 4x4 transposition in - place.
		*/
		void              transpose();

		/*
		"brief  Matrix 4x4 transpsition returns SIMDMat4x4 object.
		*/
		SIMDMat4x4        transpose() const;

		/*
		@brief  Create Identity Matrix 4x4.
		*/
		static  SIMDMat4x4        identity();

		/*
		@brief   Create upper-triangular matrix from *this, modifies  *this in - place.
		*/
		void              upper_triangular();

		/*
		"brief   Create upper-triangular matrix from *this, returns SIMDMat4x4 object
		         initialized to the content of *this.
		*/
		SIMDMat4x4             upper_triangular() const;

		/*
		@brief   Create lower-triangular matrix from *this, modifies *this in - place.
		*/
		void              lower_triangular();

		/*
		@brief   Create lower-triangualr matrix from *this, returns SIMDMat4x4 object.
		*/
		SIMDMat4x4             lower_triangular() const;

		/*
		@brief   Test for Matrix symmetricity, i.e A = AT, where T = trnspositon.
		         Returns SIMDMat4x4 filled with either one's or zeroe's or both.
		*/
		SIMDMat4x4        symmetric_matrix();

		/*
		@brief   Computes SIMDMat4x4 trace, i.e trace(SIMDMat4x4) == a00 + a11 + a22 + a33
		*/
		float              trace() const;

		/*
		@brief   Computes SIMDMat4x4 rotation in X dimension.
		*/
		static  SIMDMat4x4   x_rotation(_In_ const float);

		/*
		@brief   Computes SIMDMat4x4 rotation in Y dimension.
		*/
		static  SIMDMat4x4   y_rotation(_In_ const float);

		/*
		@brief   Computes  SIMDMat4x4 rotation in Z dimension.
		*/
		static  SIMDMat4x4   z_rotation(_In_ const float);

		/*
		@brief   Computes SIMDMat4x4 matrix scale 
		*/
		static  SIMDMat4x4   scale(_In_ const float);

		/*
		@brief   Compute SIMDMat4x4 matrix translate.
		*/
		static  SIMDMat4x4 translate(_In_ const float, _In_ const float, _In_ const float);

		/*
		@brief   Returns row 0.
		*/
		__forceinline  __m128   row_0() const;

		/*
		@brief   Returns row 1.
		*/
		__forceinline  __m128   row_1() const;

		/*
		@brief   Returns row 2.
		*/
		__forceinline  __m128   row_2() const;

		/*
		@brief   Returns row 3.
		*/
		__forceinline  __m128   row_3() const;

		/*
		@brief   Upper-Triangular Matrix Mask.
		*/
		const static __m128 UTMask[4];

		/*
		@brief   Lower-Triangular Matrix Mask.
		*/
		const static __m128 LTMask[4];

		/*
		          ***  friend operators  ***
		*/

		/*
		@brief friend operator+ returns SIMDMat4x4 object = lhs + rhs (row-wise)
		*/
		friend  SIMDMat4x4      operator+(_In_ const SIMDMat4x4 &, _In_ const SIMDMat4x4 &);

		/*
		@brief friend  operator- returns SIMDMat4x4 object = lhs - rhs (row-wise)
		*/
		friend  SIMDMat4x4      operator-(_In_ const SIMDMat4x4 &, _In_ const SIMDMat4x4 &);

		/*
		@brief  friend  operator* returns SIMDMat4x4 object = lhs * rhs (mat-mul)
		*/
		friend  SIMDMat4x4      operator*(_In_ const SIMDMat4x4 &, _In_ SIMDMat4x4 &);

		/*
		@brief  friend  operator/ returns SIMDMat4x4 object = lhs / s , where s = scalar.
		*/
		friend  SIMDMat4x4       operator/(_In_ const SIMDMat4x4 &, _In_ const float);


	private:

		/*
		@brief   Class variable  __m128  m_fMat4x4[4], represents hopefully (up to the compiler) four XMMx registers.
		*/
		__m128 m_fMat4x4[4];

	};

	/*
	     ***   Specialization for double type  ***
	*/

	template<class Args>    class  SIMDMat4x4<Args,double> {


	public:


		/*
		@brief   Default Ctor initializes m_dMat4x4[4] with zero filled entries.
		*/
		SIMDMat4x4();

		/*
		@brief   Ctor constructs an object m_dMat4x4[4] from the user passed 4x4 scalars wrapped 
		          in auxilary structure.
		*/
		SIMDMat4x4(_In_ const Args &);

		/*
		@brief    Ctor constructs an object m_dMat4x4[4] , from the user passed __m256 vector type.
		*/
		SIMDMat4x4(_In_ __m256d const &);

		/*
		@brief    Ctor constructs an object m_dMat4x4[4] , from the user passed std::initializer_list.
		*/
		SIMDMat4x4(_In_ std::initializer_list<__m256d> const &);

		/*
		@brief    Ctor constructs an object m_dMat4x4[4] , from the user passed const double.
		*/
		SIMDMat4x4(_In_ const double);

		/*
		@brief    Ctor constructs an object m_dMat4x4[4] , from the user passed reference to static type double array.
		*/
		SIMDMat4x4(_In_ double(&)[16]);

		/*
		@brief     explicit Ctor constructs an object m_dMat4x4[4] , from the user passed reference to static type float array
		*/
		explicit SIMDMat4x4(_In_ float(&)[16]);

		

		/*
		@brief    Ctor constructs an object m_dMat4x4[4] , from the user passed reference to static array of type __m256d.
		*/
		SIMDMat4x4(_In_ __m256d(&)[4]);

		/*
		@brief    Copy-Ctor.
		*/
		SIMDMat4x4(_In_ const SIMDMat4x4  &);

		/*
		@brief    Move-Ctor
		*/
		SIMDMat4x4(_In_ SIMDMat4x4 &&);

		/*
		@brief   Dtor.
		*/
		~SIMDMat4x4() = default;

	private:

		/*
		@brief   Class variable __m256d m_dMat4x4[4] , represents hopefully (up to the compiler) four YMMx registers.
		*/

		__m256d  m_dMat4x4[4];
	};

	template<typename T> struct Holder {

		T a00, a01, a02, a03;
	};

	/*
	@brief  Primary template for arguments holding.
	*/
	template<typename T> struct ArgsHolder;

	/*
	@brief float specialization.
	*/
	template<> struct  ArgsHolder<float> {

		_In_ float a00, a01, a02, a03;
		_In_ float a10, a11, a12, a13;
		_In_ float a20, a21, a22, a23;
		_In_ float a30, a31, a32, a33;
	};

	/*
	@brief  double specialization.
	*/
	template<> struct  ArgsHolder<double> {

		_In_ double a00, a01, a02, a03;
		_In_ double a10, a11, a12, a13;
		_In_ double a20, a21, a22, a23;
		_In_ double a30, a31, a32, a33;
	};

	/*
	@brief  Compile time Upper-Triangular Mask.
	*/

	template<typename T>   struct MaskUT;

	template<>  struct  MaskUT<float>{

		_In_ const __m128 v0 = _mm_set_ps(1.f,1.f,1.f,1.f);
		_In_ const __m128 v1 = _mm_set_ps(0.f, 1.f, 1.f, 1.f);
		_In_ const  __m128 v2 = _mm_set_ps(0.f, 0.f, 1.f, 1.f);
		_In_ const __m128 v3 = _mm_set_ps(0.f, 0.f, 0.f, 1.f);
	};

	template<>  struct  MaskUT<double>{

		_In_ const __m256d v0 = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
		_In_ const __m256d v1 = _mm256_set_pd(0.0, 1.0, 1.0, 1.0);
		_In_ const __m256d v2 = _mm256_set_pd(0.0, 0.0, 1.0, 1.0);
		_In_ const __m256d v3 = _mm256_set_pd(0.0, 0.0, 0.0, 1.0);
	};

	/*
	@brief  Compiler time Lower-Triangular Mask
	*/
	template<typename T>  struct MaskLT;

	/*
	     float specialization.
	*/
	template<>   struct  MaskLT<float>{

		_In_ const __m128 v0 = _mm_set_ps(1.f, 0.f, 0.f, 0.f);
		_In_ const __m128 v1 = _mm_set_ps(1.f, 1.f, 0.f, 0.f);
		_In_ const __m128 v2 = _mm_set_ps(1.f, 1.f, 1.f, 0.f);
		_In_ const __m128 v3 = _mm_set_ps(1.f, 1.f, 1.f, 1.f);
	};

	/*
	     double specialization.
	*/
	template<>    struct  MaskLT<double>{

		_In_ const __m256d v0 = _mm256_set_pd(1.0, 0.0, 0.0, 0.0);
		_In_ const __m256d v1 = _mm256_set_pd(1.0, 1.0, 0.0, 0.0);
		_In_ const __m256d v2 = _mm256_set_pd(1.0, 1.0, 1.0, 0.0);
		_In_ const __m256d v3 = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
	};

#include "GMS_sse_matrix4x4.inl"

}

}
#endif /*__GMS_SSE_MATRIX_4X4_H__*/
