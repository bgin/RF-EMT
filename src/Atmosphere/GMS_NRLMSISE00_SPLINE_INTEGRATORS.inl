
template<typename T> std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value,void>::type
atmosphere::SPLINE_INTEGRATOR<T>::splinint(_In_ const T* xa, _In_ const T* ya, _In_ const T* y2a, _In_ const int n,
_In_ const T x, _Out_ T* y)noexcept(false) {
	/*      INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
	*       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
	*       Y2A: ARRAY OF SECOND DERIVATIVES
	*       N: SIZE OF ARRAYS XA,YA,Y2A
	*       X: ABSCISSA ENDPOINT FOR INTEGRATION
	*       Y: OUTPUT VALUE
	*/
#if defined _DEBUG
	_ASSERTE(xa != nullptr && ya != nullptr && y2a != nullptr);
#else
	if (xa == nullptr || ya == nullptr || y2a == nullptr)
		throw std::runtime_error("Null pointer in: SPLINE_INTEGRATOR<T>::splinit\n");
#endif

	
	T yi(0.0);
	int klo{ 0 };
	int khi{ 1 };
	T xx, h, a, b, a2, b2;
	while ((x > xa[klo]) && (khi < n)) {
		xx = x;
		if (khi < (n - 1)) {
			if (x < xa[khi])
				xx = x;
			else
				xx = xa[khi];
		}
		h = xa[khi] - xa[klo];
		a = (xa[khi] - xx) / h;
		b = (xx - xa[klo]) / h;
		a2 = a * a;
		b2 = b * b;
		y1 += ((1.0 - a2) * ya[klo] * 0.5 + b2 * ya[khi] * 0.5 + ((-(1.0 + a2*a2)*0.25 + a2 * 0.25) *
			y2a[klo] + (b2*b2*0.25 - b2*0.5) * y2a[khi) * h * h * 0.166666666666666) * h;
		++klo;
		++khi;
	}
	*y = y1;
}

template<typename T> std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value,void>::type
atmosphere::SPLINE_INTERPOLATOR<T>::splinterp(_In_ const T* xa, _In_ const T* ya, _In_ const T* y2a,
_In_ const int n, _In_ const T, _Out_ T* y)noexcept(false) {
	/*      CALCULATE CUBIC SPLINE INTERP VALUE
	*       
	*       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
	*       Y2A: ARRAY OF SECOND DERIVATIVES
	*       N: SIZE OF ARRAYS XA,YA,Y2A
	*       X: ABSCISSA FOR INTERPOLATION
	*       Y: OUTPUT VALUE
	*/
#if defined _DEBUG
	_ASSERTE(xa != nullptr && ya != nullptr && y2a != nullptr);
#else
	if (xa == nullptr || ya == nullptr || y2a == nullptr)
		throw std::runtime_error("Null pointer in SPLINE_INTERPOLATOR<T>::splinterp\n");
#endif
	int kl0{ 0 };
	int khi{ n - 1 };
	int k;
	T h, a, b, y1;
	while ((khi - klo) > 1) {
		k = (khi + klo) / 2;
		if (xa[k] > x)
			khi = k;
		else
			klo = k;
	}
	h = xa[khi] - xa[klo];
	if (h == 0.0)
		throw std::runtime_error("Bad value in T* xa\n");
	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	y1 = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h * 0.166666666666666;
	*y = y1;
	
}

template<typename T> std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value,void>::type
atmosphere::SPLINE_CUBIC_DERIV<T>::splinderiv(_In_ const T* x, _In_ const T* y, _In_ const int n, _In_ const T yp1,
_In_ const T ypn, _Out_ T* y2)noexcept(false) {
	/*       CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
	*       
	*       X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
	*       N: SIZE OF ARRAYS X,Y
	*       YP1,YPN: SPECIFIED DERIVATIVES AT X[0] AND X[N-1]; VALUES
	*                >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
	*       Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
	*/
#if defined _DEBUG
	_ASSERTE(x != nullptr && y != nullptr && y2 != nullptr);
#else
	if (x == nullptr || y == nullptr || y2 == nullptr)
		throw std::runtime_error("Null pointer in SPLINE_CUBIC_DERIV<T>::splinderiv\n");
#endif
	T* u;
	try{
		u = new T[n];
	}
	catch (std::bad_alloc& ba){
		std::cerr << "Failed to allocated memory:" << ba.what() << std::endl;
	}
	T sig, p, qn, un;
	int i, k;
	constexpr T huge( 1.0E+30 );
	if (yp1 > huge) {
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else{
		y2[0] = -0.5;
		u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
	}
	for (i{ 1 }; i != (n - 1); ++i) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * y2[i - 1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
			(x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}
	if (ypn > huge) {
		qn = 0.0;
		un = 0.0;
	}
	else {
		qn = 0.5;
		un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2) / (x[n - 1] - x[n - 2]));
	}
	y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
	for (k = n - 2; k >= 0; --k)
		y2[k] = y2[k] * y2[k + 1] + u[k];
	delete[]u;
		
}