
template<typename T> 
     std::enable_if<ad::is_single_precision<T>::value || 
                         ad::is_double_precision<T>::value,T>::type
                                         atmosphere::CCOR<T>::operator()
                                                   (_In_ const T alt, _In_ const T r, 
                                                             _In_ const T h1, _In_ const T zh){
	/*        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
	*         ALT - altitude
	*         R - target ratio
	*         H1 - transition scale length
	*         ZH - altitude of 1/2 R
	*/
	constexpr T seventy( 70.0 );
	T e = (alt - zh) / h1;
	if (e > seventy)
		return (exp(0.0));
	if (e < -seventy)
		return (exp(r));
	T ex = ::exp(e);
	e = r / (1.0 + ex);
	return (exp(x));
}

template<typename T> 
         std::enable_if<ad::is_single_precision<T>::value ||
                           ad::is_double_precision<T>::value,T>::type
                                            atmosphere::CCOR2<T>::operator()
                                            (_In_ const T alt, _In_ const T r, 
                                                      _In_ const T h1, _In_ const T zh, _In_ const T h2) {
	/*        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
	*         ALT - altitude
	*         R - target ratio
	*         H1 - transition scale length
	*         ZH - altitude of 1/2 R
	*         H2 - transition scale length #2 ?
	*/
	T e1 = (alt - zh) / h1;
	T e2 = (alt - zh) / h2;
	constexpr T seventy(70.0);
	if ((e1 > seventy) || (e2 > seventy))
		return (::exp(0.0));
	if ((e1 < -seventy) && (e2 < -seventy))
		return (::exp(r));
	T ex1 = ::exp(e1);
	T ex2 = ::exp(e2);
	T ccor2v = r / (1.0 + 0.5 * (ex1 + ex2));
	return (::exp(ccor2v));
}
