


#include <math.h>
#include "GMS_drag_coeff.h"
#include "GMS_indices.h"

void
gms
::compute_drag_coeff(  double * __restrict Cd,
					    const int32_t nt,
					    const double * __restrict C0,
                     	const double delta0,
					    const double * __restrict Re) {
	
	for (int32_t i = 0; i != nt; ++i) {
		Cd[i] = C0[i] * (1.0 + delta0 / pow(Re[i],0.5)) * 
			            (1.0 + delta0 / pow(Re[i],0.5));
	}
}

void
gms
::compute_reynolds_number( double * __restrict Re,
						    const int32_t nt,
						    const double * __restrict Vt,
						    const double D,
						    const double * __restrict rho_f,
						    const double * __restrict d_visc,
						    const int32_t  nx,
						    const int32_t  ny,
						    const int32_t  nz) {
	
	for (int32_t i = 0; i != nt; ++i) {
		 
		for (int32_t ix = 0; ix != nx; ++ix) {
			for (int32_t iy = 0; iy != ny; ++iy) {
				for (int32_t iz = 0; iz != nz; ++iz) {
					Re[i] = (Vt[i] * D * rho_f[Ix3D(ix, nx, iy, ny, iz)]) / d_visc[Ix2D(ix,nx,iy)];
				}
			}
		}
	}
}

void
gms
::compute_davies_number( double * __restrict X,
					      const int32_t nt,
						  const double vb,
						  const double rho_b,
						  const double * __restrict rho_f,
						  const int32_t nx,
						  const int32_t ny,
						  const int32_t nz,
						  const double D,
						  const double A,
						  const double * __restrict k_visc) {

	double term1{}, term2{};
	for (int32_t i = 0; i != nt; ++i) {
		 
		for (int32_t ix = 0; ix != nx; ++ix) {
			for (int32_t iy = 0; iy != ny; ++iy) {
				for (int32_t iz = 0; iz != nz; ++iz) {
					term1 = 2.0 * vb * abs(rho_b - rho_f[Ix3D(ix,nx,iy,ny,iz)]) * 9.81 * D * D;
					term2 = A * rho_f[Ix3D(ix, nx, iy, ny, iz)] * k_visc[Ix3D(ix,nx,iy,ny,iz)] * 
							k_visc[Ix3D(ix,nx,iy,ny,iz)];
					X[i] = term1/term2;
				}
			}
		}
	}
}

void
gms
::compute_aRe( double * __restrict aRe,
			     double * __restrict bRe,
			     double * __restrict X,
			    const int32_t nt,
			    const double delta0,
			    const double * __restrict C1) {

	double term1{}, term2{};
	for (int32_t i = 0; i != nt; ++i) {
		 term1 = delta0 * delta0 * 0.25;
		 term2 = pow(pow(1.0 + C1[i] * pow(X[i], 0.5), 0.5) - 1.0, 2.0) / (pow(X[i], bRe[i]));
		 aRe[i] = term1 * term2;
	}
	
}

void
gms
::compute_bRe( double * __restrict bRe,
			    const int32_t nt,
			    const double * __restrict C1,
			    const double * __restrict X) {

	double term1{}, term2{}, term3{};
	for (int32_t i = 0; i != nt; ++i) {
		term1 = C1[i] * pow(X[i],0.5);
		term2 = 2.0 * pow(1.0 + C1[i] * pow(X[i],0.5),0.5) - 1.0;
		term3 = pow(1.0 + C1[i] * pow(X[i],0.5),0.5);
		bRe[i] = term1 / (term2 * term3);
	}
}

void
gms
::psi_interpolation( double *  __restrict psi,
					  const double * __restrict X,
					  const double * __restrict X0,
					  const int32_t nt,
					  const double k,
					  const double * __restrict Ct) {

	double term1{}, term2{};

	for (int32_t i = 0; i != nt; ++i) {

		if (X[i] < X0[i]) {
			psi[i] = 1.0;
		}
		else if (X[i] > X0[i]) {
			psi[i] = 1.0 / Ct[i];
		}
		else {
			term1 = 1.0 + pow(X[i] / X0[i],k);
			term2 = 1.0 + Ct[i] * pow(X[i] / X0[i],k);
			psi[i] = term1 / term2;
		}
		   
	}
}

void
gms
::compute_delta_bRe( double * __restrict dbRe,
					  const int32_t nt,
					  const double * __restrict X,
					  const double * __restrict X0,
					  const double * __restrict Ct,
					  const double k) {

	double term1{}, term2{}, term3{}, z{};
	for (int32_t i = 0; i != nt; ++i) {
		z = pow(X[i] / X0[i], k);
		term1 = k * Ct[i] - 1.0 * z;
		term2 = 2.0 * 1.0 + z;
		term3 = 1.0 + Ct[i] * z;
		dbRe[i] = term1 / (term2 * term3);
	}
}













