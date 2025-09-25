
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