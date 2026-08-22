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

#ifndef __GMS_ATMOSPHERE_MODEL76_H_10_03_16__
#define __GMS_ATMOSPHERE_MODEL76_H_10_03_16__


#include "GMS_AtmosphereDefs.h"
#include "GMS_AtmosPhysConstants.h"
#include "GMS_AtmMetricViscosity.h"


/*---------------------------------------------------------------------------------------------*/
/*--------------------- Adapted from  Atmosphere Model 1976 -----------------------------------*/
/*-- Based on original work of  Ralph L. Carmichael, Public Domain Aeronautical Software ------*/
/*-------      REFERENCE - U.S. Standard Atmosphere, 1976. U.S. Govt Printing Office -------   */
namespace atmosphere {

	/*
	    Local Physical Constants 
	*/

	/*template<typename T> struct PhysConsts {
		constexpr static const int NTAB{ 8 };
		constexpr static  T FT2METERS{ 0.3048 };  // mult. ft. to get meters (exact)
		constexpr static  T KELVIN2RANKINE{ 1.8 }; // mult kelvins to get deg R
		constexpr static  T PSF2NSM{ 47.880258 }; // mult lb/sq.ft to get N/sq.m
		constexpr static  T SCF2KCM{ 515.379 };  // mult slugs/cu.ft to get kg/cu.m
		constexpr static  T TZERO{ 288.15 }; // sea level temperature, kelvins
		constexpr static  T PZERO{ 101325.0 }; // sea-level pressure, N/sq.m
		constexpr static  T RHOZERO{ 1.225 };  // sea level density, kg/cu.m
		constexpr static  T AZERO{ 340.294 };   // sea-level speed of sound, m/sec
		constexpr static  T REARTH{ 6369.0 };  // Earth radius.
		constexpr static  T GMR{ 34.163195 };
		const static  T htab[NTAB];
		const static  T ttab[NTAB];
		const static  T ptab[NTAB];
		const  static T gtab[NTAB];
		
	};*/

	

	template<typename T, class PhConstants> class AtmModel76 {

	public:

		/*---------------------------------------------*/
		/*---------Constructors and Destructor---------*/
		/*---------------------------------------------*/
		AtmModel76() = delete;

		/*
		@brief  "Main" Class Ctor constructs detailed
		         atmospheric model.
		*/
		inline AtmModel76(_In_ const T);

		/*
		@brief   Copy-Ctor
		*/
		inline AtmModel76(_In_ const AtmModel76 &);

		/*
		@brief   Move-Ctor
		*/
		inline AtmModel76(_In_ AtmModel76 &&);

		/*
		@brief   Dtor = default.
		*/
		~AtmModel76()noexcept(true) = default;

		/*
		    Member and friend operators  
		*/

		/*
		 @brief      *this = rhs (copy).
		 */
		inline   auto    operator=(_In_ const AtmModel76 &)->AtmModel76 &;

		/*
		@brief        *this = rhs (move).
		*/
		inline   auto    operator=(_In_ AtmModel76 &&)->AtmModel76 &;

		/*
		@brief        subscript operator (mutable state).
		*/
		inline   auto    operator[](_In_ const int)->std::pair<std::string, T>;

		/*
		@brief        subscript operator (immutable state).
		*/
		inline   auto    operator[](_In_ const int)const->const std::pair<std::string, T>;

		/*
		@brief        operator call-forwards creation to createAtmosphere member function.
		@params       Template argument of type PhConstants
		@returns      self-reference to *this.
		*/
		inline   auto    operator()(_In_ const PhConstants &)->AtmModel76 &;

		/*
		@brief        operator<< , outputs to std::ostream.
		@params       std::ostream &
		@params       const AtmModel76Simple &
		@returns      std::ostream &
		*/
		 friend    auto    operator<<(_In_ std::ostream &, _In_ const AtmModel76 &)->std::ostream &;

		/*
		      -------    Accessor methods  --------
		*/

		/*
		@brief           Returns atmosphere data vector.
		*/
		inline   auto    getAtmData()const->std::vector<std::pair<std::string, T>>;

		/*
		@brief           Returns atmosphere read out altitude.
		*/
		inline   auto    getAltitude()const->T;


		/*
		      -------    Class   member functions   --------
		*/

		/*
		 @brief          Creates "detailed" model of Atmosphere exact up to 86km.
		*/
	    auto    createAtmosphere(_In_ const PhConstants &)->void;

		/*
		@brief           displays state of the AtmModel76 object.
		*/
		auto             displayAtmModel76() const->void;

		
	private:

		/*       
		 
		@brief    Atmospheric data table for up to 86km.
		          class member m_AtmData.
		          Data size is constant i.e 10 variables of type std::pair<std::string,T>
		*/
		std::vector<std::pair<std::string, T>> m_AtmData;
		
		/*
		@brief    Class member m_Altitude of type T scalar.
		          Altitude for which data must be computed.
		*/
		T   m_Altitude;

		
	};

	/*
	        Global function for creating and filling in static array of read-out vectors.
			This function expects pre-initialized array as an argument.
	*/

	template<typename T, class PhConstants> auto 
		createAtm76(_Inout_ atmosphere::AtmModel76<T, PhConstants>(&)[44],_In_ const PhConstants & ,_In_ const int )->void;

	/*
	         Global function for creating and filling in a vector of smart pointers pointing
			 to AtmModel76 objects.
			 This function expects AtmModel76 intialized smart pointers.
	*/
	template<typename T, class PhConstants> auto
		createAtm76(_Inout_ std::vector<std::unique_ptr<atmosphere::AtmModel76<T, PhConstants>>> &, _In_ const PhConstants &)->void;


	
#endif


























#include "GMS_AtmosphereModel76.inl"
}
#endif /*__GMS_ATMOSPHERE_MODEL76_10_03_16__*/
