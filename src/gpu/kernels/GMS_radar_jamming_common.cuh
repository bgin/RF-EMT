
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

#ifndef __GMS_RADAR_JAMMING_COMMON_CUH__
#define __GMS_RADAR_JAMMING_COMMON_CUH__

// Speed of Light
#define C 299792458.0f

// Standard temperature
#define T0 290.0f

// Effective Earth radius
#define k_e 1.333333333333333333333f

// Boltzman constant
#define k_B 1.38064852e-23f

// Earth radius
#define a_e 6378388.0f

// PI constants
#if !defined(PI)
    #define PI 3.1415926535897932384626f
#endif

#define _4PI 12.5663706143591729538506f


inline __device__
float to_dB(const float x) {
       const float tm30 = 0.000000000000000000000000000001f; //1.00000000317107685097105134714e-30
       return (10.0f*log10f(x+tm30));
}

inline __device__
float from_dB(const float x) {
      return (powf(10.0f,0.1f*x));
}

inline __device__
float integ_n_pulses(const float xtf,
                     const float xtr) {
      const float invtr = 1.0f/xtr;
      return (xtf*invtr);
}

inline __device__
float duty_cycle(const float xrho,
                 const float xtr) {
       const float invtr = 1.0f/xtr;
       return (xrho*invtr);
}

inline __device__
float radar_avg_pow(const float xPt,
                    const float Dc) {
       return (xPt*Dc);
}

inline __device__
float azimuth_bw(const float xKth,
                 const float xgamm,
		 const float xw) {
        return (xKth*(xgamm/xw));
}

inline __device__
float elevation_bw(const float xKth,
                 const float xgamm,
		 const float xh) {
     return (xKth*(xgamm/xh));
}

inline __device__
float radar_ant_gain(const float tha,
                     const float the,
		     const float xLn) {
       const float den = tha*the*xLn;
       return (_4PI/den);
}

inline __device__
float noise_density(const float xTs) {
       return (xTs*k_B);
}
































































#endif /*__GMS_RADAR_JAMMING_COMMON_CUH__*/
