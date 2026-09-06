/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#ifndef __GMS_VECTOR3D_H__
#define __GMS_VECTOR3D_H__ 161220171941



#include <iosfwd>
#include <iomanip>
#include <cstdint>

namespace file_info {



    const unsigned int gGMS_VECTOR3D_MAJOR = 1;
    const unsigned int gGMS_VECTOR3D_MINOR = 0;
    const unsigned int gGMS_VECTOR3D_MICRO = 0;
    const unsigned int gGMS_VECTOR3D_FULLVER =
      1000U*gGMS_VECTOR3D_MAJOR+100U*gGMS_VECTOR3D_MINOR+10U*gGMS_VECTOR3D_MICRO;
    const char * const pgGMS_VECTOR3D_CREATE_DATE = "16-12-2017 09:41 +00200 (SAT 16 DEC 2017 GMT+2)";
    const char * const pgGMS_VECTOR3D_BUILD_DATE  = __DATE__":"__TIME__;
    const char * const pgGMS_VECTOR3D_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_VECTOR3D_SYNOPSIS    =  "Implementation of Vector 3D (physical meaning).";
}

namespace gms {
	namespace math{

#if !defined (CHECK_VECTOR3D_LEPS)
#define CHECK_VECTOR3D_LEPS(value)  \
		do {						\
		    if ((value) <= std::numeric_limits<double>::epsilon()) \
		             (value) += std::numeric_limits<double>::epsilon();    \
		} while (0);
#endif

#if !defined (CHECK_VECTOR3D_NAN)
#define CHECK_VECTOR3D_NAN(x,y,z,bres) \
		do {					  \
		   if ((x) == std::isnan((x)) || { \
		       (y) == std::isnan((y)) ||  \
		       (z) == std::isnan((z))) {   \
			    std::cerr << "Detected: NAN at line: " << \
				 __LINE__ << " at virtual address:" << std::hex << \
				 "0x" << __FUNCTIONW__ << " in file: " << __FILE__ << "\n"; \
				 std::cerr << "x= " << (x) << \
				              "y= " << (y) << \
							  "z= " << (z) << "\n"; \
			} \
		} while (0);
#endif

#if !defined (GMS_VECTOR3D_ACCESSOR_GENERATOR)
#define GMS_VECTOR3D_ACCESSOR_GENERATOR     \
	inline double get_x() const {			\
	    return (m_x);						  \
	}										  \
	inline double get_y() const {			 \
	    return (m_y);					    \
	}										\
	inline double get_z() const {			\
		return (m_z);						\
	}										\
	inline const double * get_ptr()	 const  { \
		const double * ptr = reinterpret_cast<const double*>(&this->m_x); \
		return (ptr);														\
	}							
#endif

#if !defined (GMS_VECTOR3D_CTORS_GENERATOR)
#define GMS_VECTOR3D_CTORS_GENERATOR      \
										  \
		Vector3D();				          \
										  \
		Vector3D(_In_ const double,       \
				 _In_ const double,       \
				 _In_ const double);      \
										  \
	explicit Vector3D(_In_ const double); \
										  \
	 Vector3D(_In_ const Vector3D &);     
#endif


	class Vector3D {

			
			public:

			GMS_VECTOR3D_CTORS_GENERATOR
		    
			~Vector3D() = default;

			GMS_VECTOR3D_ACCESSOR_GENERATOR

			// Member operators

			Vector3D & operator=(const Vector3D &);

			const double operator[](const uint32_t) const;

			friend	std::ostream & operator<<(std::ostream &,
							  const Vector3D &);
			private:

			// Vector components
			double m_x;

			double m_y;

			double m_z;
		};

// Getters macro

#if !defined (V3D_X)
#define V3D_X(obj) obj.get_x()
#endif

#if !defined (V3D_Y)
#define V3D_Y(obj) obj.get_y()
#endif

#if !defined (V3D_Z)
#define V3D_Z(obj) obj.get_z()
#endif

	// Global operators

	static inline Vector3D operator+(const Vector3D &,
					 const Vector3D &);

	static inline Vector3D operator+(const Vector3D &,
					         const double);

	static inline Vector3D operator+(const double,
				         const Vector3D &);

	static inline Vector3D operator+=(Vector3D &,
				        const Vector3D &);

	static inline Vector3D operator+=(Vector3D &,
					  const double);

	static inline Vector3D operator+=(const double,
				         Vector3D &);

	static inline Vector3D operator-(const Vector3D &,
				         const Vector3D &);

	static inline Vector3D operator-(const Vector3D &,
				         const double);

	static inline Vector3D operator-(const double,
				         const Vector3D &);

	static inline Vector3D operator-=(Vector3D &,
				         const Vector3D &);

	static inline Vector3D operator-=(Vector3D &,
				         const double);

	static inline Vector3D operator-=(const double,
					   Vector3D &);

	static inline Vector3D operator*(const Vector3D &,
				         const Vector3D &);

	static inline Vector3D operator*(const Vector3D &,
				         const double);

	static inline Vector3D operator*(const double,
				         const Vector3D &);

	static inline Vector3D operator*=(Vector3D &,
				          const Vector3D &);

	static inline Vector3D operator*=(Vector3D &,
				          const double);

	static inline Vector3D operator*=(const double,
				          Vector3D &);

	static inline Vector3D operator/(const Vector3D &,
				         const Vector3D &);

	static inline Vector3D operator/(const Vector3D &,
				         const double);

	static inline Vector3D operator/(const double,
				         const Vector3D &);

	static inline Vector3D operator/=(Vector3D &,
				          const Vector3D &);

	static inline Vector3D operator/=(Vector3D &,
				          const double);

	static inline Vector3D operator/=(const double,
				          Vector3D &);

	static inline bool operator==(const Vector3D &,
			              const Vector3D &);

	static inline bool operator==(const Vector3D &,
				      const double);

	static inline bool operator==(const double,
				      const Vector3D &);

	static inline bool operator!=(const Vector3D &,
				      const Vector3D &);

	static inline bool operator!=(const Vector3D &,
				      const double);

	static inline bool operator!=(const double,
				      const Vector3D &);

	static inline bool operator>(const Vector3D &,
			             const Vector3D &);

	static inline bool operator>(const Vector3D &,
				     const double);

	static inline bool operator>(const double,
				     const Vector3D &);

	static inline bool operator<(const Vector3D &,
				     const Vector3D &);

	static inline bool operator<(const Vector3D &,
				     const double);

	static inline bool operator<(const double,
				     const Vector3D &);

	static inline bool operator>=(const Vector3D &,
				      const Vector3D &);

	static inline bool operator>=(const Vector3D &,
				      const double);

	static inline bool operator>=(const double,
				      const Vector3D &);

	static inline bool operator<=(const Vector3D &,
				      const Vector3D &);

	static inline bool operator<=(const Vector3D &,
			               const double);

	static inline bool operator<=(const double,
				      const Vector3D &);

	static inline double dot(const Vector3D &,
			         const Vector3D &);

	static inline double abs_dot(const Vector3D &,
				      const Vector3D &);

	static inline Vector3D cross(const Vector3D &,
			             const Vector3D &);

	static inline double tri_prod(const Vector3D &,
				         const Vector3D &,
				         const Vector3D &);

	static inline Vector3D dir_cos(const Vector3D &);

	static inline double norm(const Vector3D &);

	static inline Vector3D normalize(const Vector3D &);



#include "GMS_vector3D.inl"
	}
}

#endif /*__GMS_VECTOR3D_H__*/
