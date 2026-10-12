

#ifndef __GMS_THERNAL_NOISE_CUH__
#define __GMS_THERMAL_NOISE_CUH__





extern "C"
void therm_noise_range_cuda(       const float Frdr,
                                   const float Kth,
                                   const float rho,
                                   const float tf
                                   const float tf,
                                   const float tr,
				   const float th,
				   const float * __restrict d_Pt,
				   const float gamm,
				   const float * __restrict d_w,
				   const float * __restrict d_h,
				   const float Ln,
				   const float * __restrict d_Ts,
				   const float sig,
				   const float F,
				   const float Fp,
				   const float Flens,
				   const float Dx,
				   const float * __restrict d_Lt,
				   const float * __restrict d_La,
				   float * __restrict d_Rm,
				   const int nth);
				  











#endif /*__GMS_THERMAL_NOISE_CUH__*/