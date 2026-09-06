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

#ifndef GMS_TRIGFUNC_LUT_H_
#define GMS_TRIGFUNC_LUT_H_


#include "GMS_config.h"

namespace gms
{

	template<typename _Ty>  struct SinCoeffLUT
	{

	

		
		
		/*
		@brief  static const array of 8*8 Sin(2*Pi*k*i/T) elements where T=8, k=1....8, i=0...7.
		*/
		__ATTR_ALIGN__(64)
		static const   _Ty m_sSinCoeff8x8[64];

	
		/*
		@brief  static const array of 8*8 Sin(Pi*k*i/T) of half period elements where, T=8, k=1....8, i=0....7.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sSinCoeffH8x8[64];
		
		/*
		@brief  static const array of 16*16 Sin(2*Pi*k*i/T) elements where, T=16, k=1....16, i=0....15.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sSinCoeff16x16[256];

		/*
		@brief   static const array of 16*16 Sin(Pi*k*i/T) of half period elements  where, T=16, k=1....16, i=0....15.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sSinCoeffH16x16[256];

		/*
		@brief   static const array of 32*32 Sin(2*Pi*k*i/T) elements where, T=16, k=1....32, i=0....31.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sSinCoeff32x32[1024];

	};


	template<typename _Ty>  struct CosCoeffLUT
	{

	


		/*
		@brief  static constexpr array of 8*8 Cos(2*Pi*k*i/T) where T=8, k=1....8, i=0...7.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sCosCoeff8x8[64];

		/*
		@brief  static const    array  of 8*8 Cos(Pi*k*i/T) half period elements where t=8, k=1....8, i=0....7.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sCosCoeffH8x8[64];

		/*
		@brief  static const    array  of 16*16 Cos(2*Pi*k*i/T) elements where, T=16, k=1....16, i=0....15.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sCosCoeff16x16[256];

		/*
		@brief   static const array of 16*16 Cos(Pi*k*i/T) of half period elements  where, T=16, k=1....16, i=0....15.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sCosCoeffH16x16[256];

		/*
		@brief   static const array of 32*32 Sin(2*Pi*k*i/T) elements where, T=16, k=1....32, i=0....31.
		*/
		__ATTR_ALIGN__(64)
		static const  _Ty m_sCosCoeff32x32[1024];
	};

#include "GMS_trigfuncLUT.inl"
}

#endif /*GMS_TRIGFUNC_LUT_H_*/
