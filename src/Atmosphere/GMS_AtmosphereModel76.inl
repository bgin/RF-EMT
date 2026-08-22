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

/*template<typename T> const typename T atmosphere::PhysConsts<T>::htab[atmosphere::PhysConsts<T>::NTAB] = { T(0.0), T(11.0), T(20.0), T(32.0),
T(47.0), T(51.0), T(71.0), T(84.852) };

template<typename T> const typename T atmosphere::PhysConsts<T>::ttab[atmosphere::PhysConsts<T>::NTAB] = { T( 288.15 ), T( 216.65 ), T(216.65), T( 228.65),
T( 270.65 ), T( 270.65 ), T( 214.65 ), T( 186.964) };

template<typename T> const typename T atmosphere::PhysConsts<T>::ptab[atmosphere::PhysConsts<T>::NTAB] = { T(1.0), T(2.2336110E-1), T(5.4032950E-2), T(8.5666784E-3),
T(1.0945601E-3), T(6.6063531E-4), T(3.9046834E-5), T(3.68501E-6) };

template<typename T> const typename T atmosphere::PhysConsts<T>::gtab[atmosphere::PhysConsts<T>::NTAB] = { T(-6.5), T(0.0), T(1.0), T(2.8), T(0.0), T(-2.8),
T(-2.0), T(0.0) };*/


template<typename T, class PhConstants> inline atmosphere::
AtmModel76<T, PhConstants>::AtmModel76(_In_ const T alt) :
m_Altitude{ alt }
{
#if defined _DEBUG
	_ASSERT(this->m_Altitude >= static_cast<T>(-2.0) && this->m_Altitude <= static_cast<T>(46.0));
#else
	if (this->m_Altitude < static_cast<T>(-2.0) &&
		this->m_Altitude > static_cast<T>(46.0))
		throw std::runtime_error("Fatal Error in: AtmModel76<T,PhConstants>::AtmModel76(const T)");
#endif
	constexpr int const datumSize{ 10 };
	this->m_AtmData = std::vector<std::pair<std::string, T>>(datumSize);
}

template<typename T, class PhConstants> inline atmosphere::
AtmModel76<T, PhConstants>::AtmModel76(_In_ const AtmModel76 &rhs) :
m_Altitude{ rhs.m_Altitude },
m_AtmData{ rhs.m_AtmData }
{

}

template<typename T, class PhConstants> inline atmosphere::
AtmModel76<T, PhConstants>::AtmModel76(_In_  AtmModel76 &&rhs) :
m_Altitude{ std::move(rhs.m_Altitude) },
m_AtmData{ std::move(rhs.m_AtmData) }
{

}

template<typename T, class PhConstants> inline auto atmosphere::
AtmModel76<T, PhConstants>::operator=(_In_ const AtmModel76 &rhs)->atmosphere::AtmModel76<T, PhConstants> & {
	
	if (this == &rhs) return (*this);
	this->m_AtmData.operator=(rhs.m_AtmData);
	this->m_Altitude = rhs.m_Altitude;
	return (*this);
}

template<typename T, class PhConstants> inline auto atmosphere::
AtmModel76<T, PhConstants>::operator=(_In_ AtmModel76 &&rhs)->atmosphere::AtmModel76<T, PhConstants> & {

	if (this == &rhs) return (*this);
	this->m_AtmData.operator=(std::move(rhs.m_AtmData));
	this->m_Altitude = rhs.m_Altitude;
	return (*this);
}

template<typename T, class PhConstants> inline auto atmosphere::
AtmModel76<T, PhConstants>::operator()(_In_ const PhConstants &c)->atmosphere::AtmModel76<T, PhConstants> & {

	this->createAtmosphere(c);
	return (*this);
}

template<typename T, class PhConstants> inline auto atmosphere::
AtmModel76<T, PhConstants>::operator[](_In_ const int index)->std::pair<std::string, T> {
#if defined _DEBUG
	_ASSERTE(index >= 0 && index <= 9);
#else
	if(index < 0 || index > 9)
		throw std::runtime_error("Fatal Error in: AtmModel76<T,PhConstants>::operator[](int)");
#endif
	return (this->m_AtmData.operator[](index));
}

template<typename T, class PhConstants> inline auto atmosphere::
AtmModel76<T, PhConstants>::operator[](_In_ const int index)const->const std::pair<std::string, T> {
#if defined _DEBUG
	_ASSERTE(index >= 0 && index <= 9);
#else
	if(index < 0 || index > 9)
		throw std::runtime_error("Fatal Error in: AtmModel76<T,PhConstants>::operator[](int)const");
#endif
	return (this->m_AtmData.operator[](index));
}

template<typename T, class PhConstants>   auto     operator<<(_In_ std::ostream &os, _In_ const atmosphere::AtmModel76<T, PhConstants> &rhs)->std::ostream & {

	int i_size = static_cast<int>(rhs.m_AtmData.size());
	os << "Atmosphere Model-76 object state" << std::endl;
	for (int i{ 0 }; i != i_size; ++i)
		os << rhs.getAtmData()[i].first.c_str() << rhs.getAtmData()[i].second << std::endl;
	return os;
}

template<typename T, class PhConstants> inline auto      atmosphere::AtmModel76<T, PhConstants>::
getAltitude()const->T {
	return (this->m_Altitude);
}

template<typename T, class PhConstants> inline auto      atmosphere::AtmModel76<T, PhConstants>::
getAtmData()const->std::vector< std::pair<std::string, T>>{
	return (this->m_AtmData);
}


template<typename T, class PhConstants> inline auto atmosphere::
AtmModel76<T, PhConstants>::createAtmosphere(_In_ const PhConstants &c)->void {
	T delta, theta, sigma;
	constexpr static const int NTAB{ 8 };
	int i, j, k;
	T temp, pressure, density, asound;
	T viscosity, kinemat_visicosity;

	T h = this->m_Altitude * c.REARTH / (this->m_Altitude + c.REARTH);
	i = 0; j = NTAB - 1;
	do
	{
		k = (i + j) / 2;
		if (h < c.htab[k]) j = k; else i = k;

	} while (j > i + 1);

	T tgrad = c.gtab[i];
	T tbase = c.ttab[i];
	T deltah = h - c.htab[i];
	T tlocal = tbase + tgrad * deltah;
	theta = tlocal / c.ttab[0];
	if (static_cast<T>(0.0) == tgrad)
		delta = c.ptab[i] * std::exp(-c.GMR*deltah / tbase);
	else
		delta = c.ptab[i] * std::pow(tbase / tlocal, c.GMR / tgrad);
	sigma = delta / theta;
	temp = c.TZERO * theta;
	pressure = c.PZERO * delta;
	density = c.RHOZERO * sigma;
	asound = c.AZERO * std::sqrt(theta);
	atmosphere::MetricViscosity<T> mv;
	viscosity = mv.operator()(theta);
	kinemat_visicosity = viscosity / density;
	
	/* Fill in read out vector */
	this->m_AtmData.operator[](0).operator=({ std::string("Altitude"), this->m_Altitude });
	this->m_AtmData.operator[](1).operator=({ std::string("Sigma"), sigma });
	this->m_AtmData.operator[](2).operator=({ std::string("Delta"), delta });
	this->m_AtmData.operator[](3).operator=({ std::string("Theta"), theta });
	this->m_AtmData.operator[](4).operator=({ std::string("Pressure"), pressure });
	this->m_AtmData.operator[](5).operator=({ std::string("Temp"), temp });
	this->m_AtmData.operator[](6).operator=({ std::string("Density"), density });
	this->m_AtmData.operator[](7).operator=({ std::string("Sound Speed"),asound });
	this->m_AtmData.operator[](8).operator=({ std::string("Viscosity"), 1.0E+6*viscosity });
	this->m_AtmData.operator[](9).operator=({ std::string("KViscosity"), kinemat_visicosity });
}


template<typename T, class PhConstants> auto  atmosphere::AtmModel76<T, PhConstants>::
displayAtmModel76()const->void {
	std::printf("****AtmModel76 object state****:\n\n");
	std::printf("units of alt are Km; pressure N/sq.m.; density is kg/cu.m.\n");
	std::printf("speed of sound is m/s; viscosity is kg per m-s\n");
	std::printf("%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n", "alt", "sigma", "delta", "theta", "press", "temp", "dens", "v", "visc", "k.visc");

	std::printf("%-10.3f%-10.3f%-10.3f%-10.3f%-10.3f%-10.3f%-10.3f%-10.3f%-10.3f%-10.3f\n", this->m_AtmData[0].second, this->m_AtmData[1].second,
		this->m_AtmData[2].second, this->m_AtmData[3].second, this->m_AtmData[4].second, this->m_AtmData[5].second, this->m_AtmData[6].second,
		this->m_AtmData[7].second, this->m_AtmData[8].second, this->m_AtmData[9].second);
	
}




template<typename T, class PhConstants> auto 
createAtm76(_Inout_ atmosphere::AtmModel76<T, PhConstants>(&Sys)[44], _In_ const PhConstants &c, _In_ const int datumSize)->void {
	constexpr int SysLength{ 44U };
#if defined _DEBUG
	_ASSERTE(datumSize == SysLength);
#else
	if(datumSize != SysLength) throw std::runtime_error("Fatal Error in createAtm76(AtmModel76(&)[44])");
#endif
	// Static array Sys must be initialized by the caller.
	for (int i{ 0 }; i != datumSize; ++i)
		Sys[i].createAtmosphere(c);
}


 template<typename T, class PhConstants> auto
	 createAtm76(_Inout_ std::vector<std::unique_ptr<atmosphere::AtmModel76<T, PhConstants>>> &vAtmSys, _In_ const PhConstants &c)->void {
		 constexpr std::size_t SysSize{ 44U };
#if defined _DEBUG
		 _ASSERTE(vAtmSys.size() == SysSize);
#else
		 if(vAtmSys.size() != SysSize) throw std::runtime_error("Fatal Error in createAtm76(std::vector<...>)");
#endif

		 // vector must contain initialized smart pointers.
		 for (std::size_t i{ 0 }; i != vAtmSys.size(); ++i)
			 vAtmSys.operator[](i).operator->()->createAtmosphere(c);
	 }
	 
