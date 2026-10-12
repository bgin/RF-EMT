

#include <cuda_runtime.h>
#include <cuda.h>
#include "GMS_thermal_noise.cuh"
#include "GMS_jamming_common.cuh"


__global__ void
therm_noise_range_kernel(const float Frdr,
                           const float Kth,
                           const float rho,
                           const float tf,
                           const float tr,
	                   const float th,
		           const float * __restrict Pt,
		           const float gamm,
			   const float * __restrict w,
			   const float * __restrict h,
			   const float Ln,
			   const float * __restrict Ts,
			   const float sig,
			   const float F,
			   const float Fp,
		           const float Flens,
			   const float Dx,
			   const float * __restrict Lt,
			   const float * __restrict La,
			   float * __restrict Rm,
			   const int nth) {

     int tid = blockDim.x*blockIdx.x+threadIdx.x;
     if(tid < nth) {
        const float dc  = duty_cycle(rho,tr);
        const float Pav = radar_avg_pow(Pt[tid],dc);
        const float ag  = radar_ant_gain(azimuth_bw(Kth,gamm,
                                                    h[tid]),
                                         elevation_bw(Kth,gamm,
                                                    w[tid]),Ln);
        const float N0  = noise_density(Ts[tid]);
        const float den = 1984.4017075391884912304967f*N0*Dx*Lt[tid]*La[tid];
        const float t1  = gamm*gamm;
        const float t2  = Pav*tf;
        const float t3  = ag*ag;
        const float t4  = sig*Frdr*Frdr*Fp*Fp;
        const float t5  = F*F*F*F*Flens*Flens;
        const float num = t1*t2*t3*t4*t5;
        const float rat = num/den;
        Rm[tid]         = powf(rat,0.25f);
     }
}


 therm_noise_range_cuda(           const float Frdr,
                                   const float Kth,
                                   const float rho,
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
				   const int nth) {
				  
          int threadsBlock = 256;
          int blocksGrid  = (nth + threadsBlock - 1) / threadsBlock;
          therm_noise_range_kernel<<<blockGrid,threadsBlock>>>(Frdr,
                                                                 Kth,
                                                                 rho,
                                                                 tf,
                                                                 tr,
                                                                 th,
                                                                 d_Pt,
                                                                 gamm,
                                                                 d_w,
                                                                 d_h,
                                                                 Ln,
                                                                 d_ts,
                                                                 sig,
                                                                 F,
                                                                 Fp,
                                                                 Flens,
                                                                 Dx,
                                                                 d_Lt,
                                                                 d_La,
                                                                 d_Rm,
                                                                 nth);


}
