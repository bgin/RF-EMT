
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

template<typename T> std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value,T>::type
atmosphere::SCALH<T>::operator()(_In_ const T alt, _In_ const T xm, _In_ const T temp) {

	T rgas = 831.4;
	T g = ang::gsurf / (1.0 + alt / ang::re) * (1.0 + alt / ang::re);
	g = rgas * temp / (g * xm);
	return g;
}

template<typename T> std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value,T>::type
atmosphere::DNET<T>::operator()(_In_ T dd, _In_ const T dm, _In_ const T zhm, _In_ const T xmm, _In_ const T xm){

	T a = zhm / (xmm - xm);
	if (!((dm > 0.0) && (dd > 0.0))) {
		std::cerr << "DNET log args error:" << dm << dd << xm << std::endl;
		if ((dd == 0.0) && (dm == 0.0))
			dd = 1.0;
		if (dm == 0.0)
			return dd;
		if (dd == 0.0)
			return dm;
	}
	T ylog = a * log(dm / dd);
	if (ylog < -10.0)
		return dd;
	if (ylog > 10.0)
		return dm;
	a = dd * pow((1.0 + exp(ylog)), (1.0 / a));
	return a;
}
