template<typename T, class SPLINE , class SPLINT , class SPLINI>
	auto atmosphere::DENSM<T, SPLINE, SPLINT, SPLINI>::zeta(_In_ const T zz, _In_ const T z1)->T {

		return ((zz - z1) * (atmosphere::nrlmsise00_globals::re + z1) / (atmosphere::nrlmsise00_globals::re + zz));
	}

template<typename T, class SPLINE , class SPLINT  
	, class SPLINI >
	std::enable_if<ad::is_single_precision<T>::value || ad::is_double_precision<T>::value,T>::type
	atmosphere::DENSM<T, SPLINE, SPLINT, SPLINI>::densprof(_In_ const T alt, _In_ const T d0, _In_ const T xm, _In_ T* tz,
	_In_ const int mn3, _In_ const T* zn3, _In_ const T* tn3, _In_ const T* tgn3, _In_ const int mn2, _In_ const T* zn2,
	_In_ const T* tn2, _In_ const T* tgn2) {

		T xs[10], ys[10], y2out[10];
		constexpr T rgas = 831.4;
		constexpr int n_threads{ 4 };
		T z, z1, z2, t1, t2, zg, zgdif;
		T yd1, yd2;
		T x, y, yi;
		T expl, gamm, glb;
		T densm_temp = d0;
		int mn, k;
		if (alt > zn2[0]) {
			if (xm == 0.0)
				return *tz;
			else
				return d0;
		}
		/* STRATOSPHERE/MESOSPHERE TEMPERATURE */
		if (alt > zn2[mn2 - 1])
			z = alt;
		else
			z = zn2[mn2 - 1];
		mn = mn2;
		z1 = zn2[0];
		z2 = zn2[mn - 1];
		t1 = tn[0];
		t2 = tn2[mn - 1];
		zg = zeta(z, z1);
		zgdif = zeta(z2, z1);
		/* set up spline nodes */
#if defined USE_OPENMP
		omp_set_num_threads(n_threads);
#pragma omp parallel for private(k) default(shared)
		for (k{ 0 }; k < mn; ++k) {
			xs[k] = zeta(zn2[k], z1) / zgdif;
			ys[k] = 1.0 / tn2[k];
		}
#else
		for (k{ 0 }; k < mn; ++k) {
			xs[k] = zeta(zn2[k], z1) / zgdif;
			ys[k] = 1.0 / tn2[k];
		}
#endif
		T term1 = (atmosphere::nrlmsise00_globals::re + z2) / (atmosphere::nrlmsise00_globals::re + z1);
		yd1 = -tgn2[0] / (t1 * t1) * zgdif;
		yd2 = -tgn2[1] / (t2 * t2) * zgdif * term1 * term1;
		/* calculate spline coefficients */
		SPLINE::splinderiv(xs, ys, mn, yd1, yd2, y2out);
		x = zg / zgdif;
		SPLINT::splinterp(xs, ys, y2out, mn, x, &y);

		/* temperature at altitude*/
		*tz = 1.0 / y;
		if (xm != 0.0) {
			/* stratosphere/mesosphere density */
			glb = atmosphere::nrlmsise00_globals::gsurf / ((1.0 + z1 / atmosphere::nrlmsise00_globals::re) * (1.0 + z1 / atmosphere::nrlmsise00_globals::re));
			gamm = xm * glb * zgdif / rgas;
			/*
			Integrate  temperature profile */
			SPLINI::splinint(xs, ys, y2out, mn, x, &yi);
			expl = gamm * yi;
			if (expl > 50.0)
				expl = 50.0;
			/* Density at altitude */
			densm_temp = densm_temp * (t1 / *tz) * exp(-expl);
		}
		if (alt > zn3[0]) {
			if (xm == 0.0)
				return *tz;
			else
				return densm_temp;
		}
		/* troposphere/stratosphere temperature */
		z = alt;
		mn = mn3;
		z1 = zn3[0];
		z2 = zn3[mn - 1];
		t1 = tn3[0];
		t2 = tn3[mn - 1];
		zg = zeta(z, z1);
		zgdif = zeta(z2, z1);
		/* set up spline nodes */
#if defined USE_OPENMP
		omp_set_num_threads(n_threads);
#pragma omp parallel for private(k) default(shared)
		for (k{ 0 }; k < mn; ++k) {
			xs[k] = zeta(zn3[k], z1) / zgdif;
			ys[k] = 1.0 / tn3[k];
		}
#else
		for (k{ 0 }; k < mn; ++k) {
			xs[k] = zeta(zn3[k], z1) / zgdif;
			ys[k] = 1.0 / tn3[k];
		}
#endif
		T term2 = (atmosphere::nrlmsise00_globals::re + z2) / (atmosphere::nrlmsise00_globals::re + z1);
		yd1 =- tgn3[0] / (t1 * t1) * zgdif;
		yd2 = -tgn3[1] / (t2 * t2) * zgdif * term2 * term2;

		/* calculate spline coeffs */
		SPLINE::splinderiv(xs, ys, mn, yd1, yd2, y2out);
		x = zg / zgdif;
		SPLINT::splinterp(xs, ys, y2out, mn, x, &y);

		/* temperature at altitude */
		*tz = 1.0 / y;
		if (xm != 0.0) {
			/* tropospheric/stratospheric density */
			glb = atmosphere::nrlmsise00_globals::gsurf / ((1.0 + z1 / atmosphere::nrlmsise00_globals::re) * (
				1.0 + z1 / atmosphere::nrlmsise00_globals::re));
			gamm = xm * glb * zgdif / rgas;
			/* Integrate temperature profile */
			SPLINI::splinint(xs, ys, y2out, mn, x, &yi);
			expl = gamm * yi;
			if (expl > 50.0)
				expl = 50.0;

			/* Density at altitude */
			densm_temp = densm_temp * (t1 / *tz) * exp(-expl);
		}
		if (xm == 0.0)
			return *tz;
		return densm_temp;
			
	}

	template<typename T, class SPLINE, class SPLINT, class SPLINI> std::enable_if<ad::is_single_precision<T>::value || 
		ad::is_double_precision<T>::value, T>::type atmosphere::DENSM<T, SPLINE, SPLINT, SPLINI>::densuprof(_In_ const T alt, _In_ const T dlb, _In_ const T tinf,
		_In_ const T tlb, _In_ const T xm, _In_ const T alpha, _In_ const T* tz, _In_ const T zlb, _In_ const T s2,
		_In_ const int mn1, _In_ const T* zn1, _In_ const T* tn1, _In_ const T* tgn1) {

			T yd2, yd1, x, y;
			constexpr T rgas = 831.4;
			T densu_temp = 1.0;
			T za, z, zg2, tt, ta;
			T dta, z1, z2, t1, t2, zg, zgdif;
			int mn, k;
			const int nThreads{ 4 };
			T glb, expl, y1, densa, gamma, gamm;
			T xs[5], ys[5], y2out[5];
			/* joining altitudes of Bates and spline */
			za = zn1[0];
			if (alt > za)
				z = alt;
			z = za;
			/* geopotential altitude difference from ZLB */
			zg2 = zeta(z, zlb);
			/* Bates Temperature*/
			tt = tinf - (tinf - tlb) * exp(-s2 * zg2);
			ta = tt;
			*tz = tt;
			densu_temp = *tz;

			if (alt < za) {
				/* calculate temperature below ZA
				* temperature gradient at ZA from Bates profile */
				dta = (tinf - ta) * s2 * (((atmosphere::nrlmsise00_globals::re + zlb) / (atmosphere::nrlmsise00_globals::re + za)) *
					((atmosphere::nrlmsise00_globals::re + zlb) / (atmosphere::nrlmsise00_globals::re + za)));
				tgn1[0] = dta;
				tn1[0] = ta;
				if (alt > zn1[mn1 - 1])
					z = alt;
				z = zn1[mn1 - 1];
				mn = mn1;
				z1 = zn1[0];
				z2 = zn1[mn - 1];
				t1 = tn1[0];
				t2 = tn1[mn - 1];
				/* geopotential difference from z1 */
				zg = zeta(z, z1);
				zgdif = zeta(z2, z1);
				T inv_zgdif(1.0 / zgdif);
#if defined USE_OPENMP
				omp_set_num_threads(nThreads);
#pragma omp parallel for private(k) default(shared)
				for (k{ 0 }; k < mn; ++k) {
					xs[k] = zeta(zn1[k], z1) * inv_zgdif;
					ys[k] = 1.0 / tn1[k];
				}
#else
				for (k{ 0 }; k < mn; ++k) {
					xs[k] = zeta(zn1[k], z1) * inv_zgdif;
					ys[k] = 1.0 / tn1[k];
				}
#endif
				/* end node derivatives */
				yd1 = -tgn1[0] / (t1 * t1) * zgdif;
				yd2 = -tgn1[1] / (t2 * t2) * zgdif * (((atmosphere::nrlmsise00_globals::re + z2) / (atmosphere::nrlmsise00_globals::re + z1)) *
					((atmosphere::nrlmsise00_globals::re + z2) / (atmosphere::nrlmsise00_globals::re + z1)));
				/* calculate spline coefficients */
				SPLINE::splinderiv(xs, ys, mn, yd1, yd2, y2out);
				x = zg / zgdif;
				SPLINT::splinterp(xs, ys, y2out, mn, x, &y);
				/* temperature at altitude */
				*tz = 1.0 / y;
				densu_temp = *tz;
			}
			if (xm == 0.0)
				return densu_temp;

			/* calculate denisty above za */
			glb = atmosphere::nrlmsise00_globals::gsurf * ((1.0 + zlb / atmosphere::nrlmsise00_globals::re) * (1.0 + zlb / atmosphere::nrlmsise00_globals::re));
			gamma = xm * glb / (s2 * rgas * tinf);
			expl = exp(-s2 * gamma * zg2);
			if (expl > 50.0)
				expl = 50.0;
			if (tt <= 0.0)
				expl = 50.0;
			/* density at altitude */
			densa = dlb * pow((tlb / tt), (1.0 + alpha + gamma)) * expl;
			densu_temp = densa;
			if (alt >= za) /* potential error !!*/
				return densu_temp;

			/* calculate density below za */
			glb = atmosphere::nrlmsise00_globals::gsurf / ((1.0 + z1 / atmosphere::nrlmsise00_globals::re) * (1.0 + z1 /
				atmosphere::nrlmsise00_globals::re));
			gamm = xm * glb * zgdif / rgas;
			/* Integrate spline temperatures */
			SPLINI::splinint(xs, ys, y2out, mn, x, &yi);
			expl = gamm * yi;
			if (expl > 50.0)
				expl = 50.0;
			if (*tz <= 0.0)
				expl = 50.0

				/* density at altitude */
				densu_temp = densu_temp * pow((t1 / *tz), (1.0 + alpha)) * exp(-expl);
			return densu_temp;
		}