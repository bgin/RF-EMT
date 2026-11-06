
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


#ifndef __GMS_MODULE_BL_MYNN_C_IFACE_H__
#define __GMS_MODULE_BL_MYNN_C_IFACE_H__

#include <cstdint>

namespace file_info 
{

     const unsigned int GMS_MODULE_BL_MYNN_C_IFACE_MAJOR = 1;
     const unsigned int GMS_MODULE_BL_MYNN_C_IFACE_MINOR = 0;
     const unsigned int GMS_MODULE_BL_MYNN_C_IFACE_MICRO = 0;
     const unsigned int GMS_MODULE_BL_MYNN_C_IFACE_FULLVER =
       1000U*GMS_MODULE_BL_MYNN_C_IFACE_MAJOR+100U*GMS_MODULE_BL_MYNN_C_IFACE_MINOR+
       10U*GMS_MODULE_BL_MYNN_C_IFACE_MICRO;
     const char * const GMS_MODULE_BL_MYNN_C_IFACE_CREATION_DATE = "05-11-2026 10:03 AM +00200 (THR 11 NOV 2026 GMT+2)";
     const char * const GMS_MODULE_BL_MYNN_C_IFACE_BUILD_DATE    = __DATE__ " " __TIME__;
     static const char  GMS_MODULE_BL_MYNN_C_IFACE_AUTHOR[]      = "Programmer: Bernard Gingold, beniekg@gmail.com"; 
     const char * const GMS_MODULE_BL_MYNN_C_IFACE_SYNOPSIS      = "C-Interface to Fortran 90 Module module_bl_mynn.";

}

/**********************************************
C-Interface to Fortran 90 Module module_bl_mynn.
***********************************************/

#if defined(__cplusplus)

extern "C" {

/**********************************************************
    Module module_bl_mynn variables and its subroutines
   declarations.
***********************************************************/

/*****************************************
 *      Module parameters declarations.  *
 *****************************************/

// ! The parameters below depend on stability functions of module_sf_mynn.
float MODULE_BL_MYNN_mp_cphm_st;

float MODULE_BL_MYNN_mp_cphm_unst;

float MODULE_BL_MYNN_mp_cphh_st;

float MODULE_BL_MYNN_mp_cphh_unst;

float MODULE_BL_MYNN_mp_xlvcp;

float MODULE_BL_MYNN_mp_xlscp;

float MODULE_BL_MYNN_mp_ev;

float MODULE_BL_MYNN_mp_rd;

float MODULE_BL_MYNN_mp_rk;

float MODULE_BL_MYNN_mp_svp11;

float MODULE_BL_MYNN_mp_p608;

float MODULE_BL_MYNN_mp_ep_3;

// reference temperature (K)
float MODULE_BL_MYNN_mp_tref;

// for total water conversion, Tripoli and Cotton (1981).
float MODULE_BL_MYNN_mp_TKmin;

float MODULE_BL_MYNN_mp_tv0;

float MODULE_BL_MYNN_mp_tv1;

float MODULE_BL_MYNN_mp_gtr;

// Closure constants.
float MODULE_BL_MYNN_mp_vk;

float MODULE_BL_MYNN_mp_pr;

float MODULE_BL_MYNN_mp_g1;

float MODULE_BL_MYNN_mp_b1;

float MODULE_BL_MYNN_mp_b2;

float MODULE_BL_MYNN_mp_c2;

float MODULE_BL_MYNN_mp_c3;

float MODULE_BL_MYNN_mp_c4;

float MODULE_BL_MYNN_mp_c5;

float MODULE_BL_MYNN_mp_a1;

float MODULE_BL_MYNN_mp_c1;

float MODULE_BL_MYNN_mp_a2;

float MODULE_BL_MYNN_mp_g2;

float MODULE_BL_MYNN_mp_cc2;

float MODULE_BL_MYNN_mp_cc3;

float MODULE_BL_MYNN_mp_e1c;

float MODULE_BL_MYNN_mp_e2c;

float MODULE_BL_MYNN_mp_e3c;

float MODULE_BL_MYNN_mp_e4c;

float MODULE_BL_MYNN_mp_e5c;

/*
        ! Constants for min tke in elt integration (qmin), max z/L in els
(zmax), ! and factor for eddy viscosity for TKE (Kq = Sqfac*Km):
*/
float MODULE_BL_MYNN_mp_qmin;

float MODULE_BL_MYNN_mp_zmax;

float MODULE_BL_MYNN_mp_Sqfac;

/*
        ! Note that the following mixing-length constants are now specified in
   mym_length !      &cns=3.5, alp1=0.23, alp2=0.3, alp3=3.0, alp4=10.0,
   alp5=0.4
*/

/*
        ! Constants for gravitational settling
*/
float MODULE_BL_MYNN_mp_gno;

float MODULE_BL_MYNN_mp_gpw;

float MODULE_BL_MYNN_mp_qcgmin;

float MODULE_BL_MYNN_mp_qkemin;

/*
         Constants for cloud PDF (mym_condensation)
*/
float MODULE_BL_MYNN_mp_rr2;

float MODULE_BL_MYNN_mp_rrp;

/*
        'parameters' for Poisson distribution (StEM EDMF scheme)
*/
float MODULE_BL_MYNN_mp_zero;

float MODULE_BL_MYNN_mp_half;

float MODULE_BL_MYNN_mp_one;

float MODULE_BL_MYNN_mp_two;

/*
        !Use Canuto/Kitamura mod (remove Ric and negative TKE) (1:yes, 0:no)
!For more info, see Canuto et al. (2008 JAS) and Kitamura (Journal of the
!Meteorological Society of Japan, Vol. 88, No. 5, pp. 857-864, 2010).
!Note that this change required further modification of other parameters
!above (c2, c3). If you want to remove this option, set c2 and c3 constants
!(above) back to NN2009 values (see commented out lines next to the
!parameters above). This only removes the negative TKE problem
!but does not necessarily improve performance - neutral impact.
*/
float MODULE_BL_MYNN_mp_CKmod;

/*
        !Use Ito et al. (2015, BLM) scale-aware (0: no, 1: yes). Note that this
also has impacts !on the cloud PDF and mass-flux scheme, using Honnert et al.
(2011) similarity function !for TKE in the upper PBL/cloud layer.
*/
float MODULE_BL_MYNN_mp_scaleaware;

/*
   Temporary switch to deactivate the mixing of chemical species (already done
   when WRF_CHEM = 1)
*/
int MODULE_BL_MYNN_mp_bl_mynn_mixchem;

/*
        Constants used for empirical calculations of saturation
! vapor pressures (in function "esat") and saturation mixing ratios
! (in function "qsat"), reproduced from module_mp_thompson.F,
! v3.6
*/
float MODULE_BL_MYNN_mp_J0;

float MODULE_BL_MYNN_mp_J1;

float MODULE_BL_MYNN_mp_J2;

float MODULE_BL_MYNN_mp_J3;

float MODULE_BL_MYNN_mp_J4;

float MODULE_BL_MYNN_mp_J5;

float MODULE_BL_MYNN_mp_J6;

float MODULE_BL_MYNN_mp_J7;

float MODULE_BL_MYNN_mp_J8;

float MODULE_BL_MYNN_mp_K0;

float MODULE_BL_MYNN_mp_K1;

float MODULE_BL_MYNN_mp_K2;

float MODULE_BL_MYNN_mp_K3;

float MODULE_BL_MYNN_mp_K4;

float MODULE_BL_MYNN_mp_K5;

float MODULE_BL_MYNN_mp_K6;

float MODULE_BL_MYNN_mp_K7;

float MODULE_BL_MYNN_mp_K8;

int MODULE_BL_MYNN_mp_mynn_level;

unsigned char MODULE_BL_MYNN_mp_mynn_meassage[128];

int MODULE_BL_MYNN_mp_kdebug;

/*
! **********************************************************************
! *   An improved Mellor-Yamada turbulence closure model               *
! *                                                                    *
! *                                   Aug/2005  M. Nakanishi (N.D.A)   *
! *                        Modified:  Dec/2005  M. Nakanishi (N.D.A)   *
! *                                             naka@nda.ac.jp         *
! *                                                                    *
! *   Contents:                                                        *
! *     1. mym_initialize  (to be called once initially)               *
! *        gives the closure constants and initializes the turbulent   *
! *        quantities.                                                 *
! *    (2) mym_level2      (called in the other subroutines)           *
! *        calculates the stability functions at Level 2.              *
! *    (3) mym_length      (called in the other subroutines)           *
! *        calculates the master length scale.                         *
! *     4. mym_turbulence                                              *
! *        calculates the vertical diffusivity coefficients and the    *
! *        production terms for the turbulent quantities.              *
! *     5. mym_predict                                                 *
! *        predicts the turbulent quantities at the next step.         *
! *     6. mym_condensation                                            *
! *        determines the liquid water content and the cloud fraction  *
! *        diagnostically.                                             *
! *                                                                    *
! *             call mym_initialize                                    *
! *                  |                                                 *
! *                  |<----------------+                               *
! *                  |                 |                               *
! *             call mym_condensation  |                               *
! *             call mym_turbulence    |                               *
! *             call mym_predict       |                               *
! *                  |                 |                               *
! *                  |-----------------+                               *
! *                  |                                                 *
! *                 end                                                *
! *                                                                    *
! *   Variables worthy of special mention:                             *
! *     tref   : Reference temperature                                 *
! *     thl     : Liquid water potential temperature               *
! *     qw     : Total water (water vapor+liquid water) content        *
! *     ql     : Liquid water content                                  *
! *     vt, vq : Functions for computing the buoyancy flux             *
! *                                                                    *
! *     If the water contents are unnecessary, e.g., in the case of    *
! *     ocean models, thl is the potential temperature and qw, ql, vt   *
! *     and vq are all zero.                                           *
! *                                                                    *
! *   Grid arrangement:                                                *
! *             k+1 +---------+                                        *
! *                 |         |     i = 1 - nx                         *
! *             (k) |    *    |     j = 1 - ny                         *
! *                 |         |     k = 1 - nz                         *
! *              k  +---------+                                        *
! *                 i   (i)  i+1                                       *
! *                                                                    *
! *     All the predicted variables are defined at the center (*) of   *
! *     the grid boxes. The diffusivity coefficients are, however,     *
! *     defined on the walls of the grid boxes.                        *
! *     # Upper boundary values are given at k=nz.                   *
! *                                                                    *
! *   References:                                                      *
! *     1. Nakanishi, M., 2001:                                        *
! *        Boundary-Layer Meteor., 99, 349-378.                        *
! *     2. Nakanishi, M. and H. Niino, 2004:                           *
! *        Boundary-Layer Meteor., 112, 1-31.                          *
! *     3. Nakanishi, M. and H. Niino, 2006:                           *
! *        Boundary-Layer Meteor., (in press).                         *
! *     4. Nakanishi, M. and H. Niino, 2009:                           *
! *        Jour. Meteor. Soc. Japan, 87, 895-912.                      *
! **********************************************************************
*/

/*
        !     SUBROUTINE  mym_initialize:
!
!     Input variables:
!       iniflag         : <>0; turbulent quantities will be initialized
!                         = 0; turbulent quantities have been already
!                              given, i.e., they will not be initialized
!       nx, ny, nz      : Dimension sizes of the
!                         x, y and z directions, respectively
!       tref            : Reference temperature                      (K)
!       dz(nz)        : Vertical grid spacings                     (m)
!                         # dz(nz)=dz(nz-1)
!       zw(nz+1)        : Heights of the walls of the grid boxes     (m)
!                         # zw(1)=0.0 and zw(k)=zw(k-1)+dz(k-1)
!       h(nx,ny)        : G^(1/2) in the terrain-following coordinate
!                         # h=1-zg/zt, where zg is the height of the
!                           terrain and zt the top of the model domain
!       pi0(nx,my,nz) : Exner function at zw*h+zg             (J/kg K)
!                         defined by c_p*( p_basic/1000hPa )^kappa
!                         This is usually computed by integrating
!                         d(pi0)/dz = -h*g/tref.
!       rmo(nx,ny)      : Inverse of the Obukhov length         (m^(-1))
!       flt, flq(nx,ny) : Turbulent fluxes of sensible and latent heat,
!                         respectively, e.g., flt=-u_*Theta_* (K m/s)
!! flt - liquid water potential temperature surface flux
!! flq - total water flux surface flux
!       ust(nx,ny)      : Friction velocity                        (m/s)
!       pmz(nx,ny)      : phi_m-zeta at z1*h+z0, where z1 (=0.5*dz(1))
!                         is the first grid point above the surafce, z0
!                         the roughness length and zeta=(z1*h+z0)*rmo
!       phh(nx,ny)      : phi_h at z1*h+z0
!       u, v(nx,nz,ny): Components of the horizontal wind        (m/s)
!       thl(nx,nz,ny)  : Liquid water potential temperature
!                                                                    (K)
!       qw(nx,nz,ny)  : Total water content Q_w                (kg/kg)
!
!     Output variables:
!       ql(nx,nz,ny)  : Liquid water content                   (kg/kg)
!       v?(nx,nz,ny)  : Functions for computing the buoyancy flux
!       qke(nx,nz,ny) : Twice the turbulent kinetic energy q^2
!                                                              (m^2/s^2)
!       tsq(nx,nz,ny) : Variance of Theta_l                      (K^2)
!       qsq(nx,nz,ny) : Variance of Q_w
!       cov(nx,nz,ny) : Covariance of Theta_l and Q_w              (K)
!       el(nx,nz,ny)  : Master length scale L                      (m)
!                         defined on the walls of the grid boxes
!
!     Work arrays:        see subroutine mym_level2
!       pd?(nx,nz,ny) : Half of the production terms at Level 2
!                         defined on the walls of the grid boxes
!       qkw(nx,nz,ny) : q on the walls of the grid boxes         (m/s)
!
!     # As to dtl, ...gh, see subroutine mym_turbulence.
!
!-------------------------------------------------------------------
*/

void MODULE_BL_MYNN_mp_MYM_INITIALIZE(int *, int *, float *, float *, float *,
                                      float *, float *, float *, float *,
                                      float *, float *, float *, float *,
                                      float *, float *, float *, float *,
                                      float *, float *, float *, int *, float *,
                                      float *, int *);

/*
       !
! ==================================================================
!     SUBROUTINE  mym_level2:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:
!       dtl(nx,nz,ny) : Vertical gradient of Theta_l             (K/m)
!       dqw(nx,nz,ny) : Vertical gradient of Q_w
!       dtv(nx,nz,ny) : Vertical gradient of Theta_V             (K/m)
!       gm (nx,nz,ny) : G_M divided by L^2/q^2                (s^(-2))
!       gh (nx,nz,ny) : G_H divided by L^2/q^2                (s^(-2))
!       sm (nx,nz,ny) : Stability function for momentum, at Level 2
!       sh (nx,nz,ny) : Stability function for heat, at Level 2
!
!       These are defined on the walls of the grid boxes.
!
*/

void MODULE_BL_MYNN_mp_MYM_LEVEL2(int *, int *, float *, float *, float *,
                                  float *, float *, float *, float *, float *,
                                  float *, float *, float *, float *, float *,
                                  float *, float *);

/*
! ==================================================================
!     SUBROUTINE  mym_length:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:   see subroutine mym_initialize
!
!     Work arrays:
!       elt(nx,ny)      : Length scale depending on the PBL depth    (m)
!       vsc(nx,ny)      : Velocity scale q_c                       (m/s)
!                         at first, used for computing elt
!
!     NOTE: the mixing lengths are meant to be calculated at the full-
!           sigmal levels (or interfaces beween the model layers).
!
*/

void MODULE_BL_MYNN_mp_MYM_LENGTH(int *, int *, float *, float *, float *,
                                  float *, float *, float *, float *, float *,
                                  float *, float *, float *, float *, float *,
                                  float *, float *, int *, float *, float *,
                                  int *);

/*
! ==================================================================
SUBROUTINE boulac_length0(k,kts,kte,zw,dz,qtke,theta,lb1,lb2)
!
!    NOTE: This subroutine was taken from the BouLac scheme in WRF-ARW
!          and modified for integration into the MYNN PBL scheme.
!          WHILE loops were added to reduce the computational expense.
!          This subroutine computes the length scales up and down
!          and then computes the min, average of the up/down
!          length scales, and also considers the distance to the
!          surface.
!
!      dlu = the distance a parcel can be lifted upwards give a finite
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------
*/

void MODULE_BL_MYNN_mp_BOULAC_LENGTH0(int *, int *, int *, float *, float *,
                                      float *, float *, float *, float *);

/*
! ==================================================================
SUBROUTINE boulac_length(kts,kte,zw,dz,qtke,theta,lb1,lb2)
!
!    NOTE: This subroutine was taken from the BouLac scheme in WRF-ARW
!          and modified for integration into the MYNN PBL scheme.
!          WHILE loops were added to reduce the computational expense.
!          This subroutine computes the length scales up and down
!          and then computes the min, average of the up/down
!          length scales, and also considers the distance to the
!          surface.
!
!      dlu = the distance a parcel can be lifted upwards give a finite
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------
*/

void MODULE_BL_MYNN_mp_BOULAC_LENGTH(int *, int *, float *, float *, float *,
                                     float *, float *, float *);

/*
! ==================================================================
!     SUBROUTINE  mym_turbulence:
!
!     Input variables:    see subroutine mym_initialize
!       levflag         : <>3;  Level 2.5
!                         = 3;  Level 3
!
!     # ql, vt, vq, qke, tsq, qsq and cov are changed to input variables.
!
!     Output variables:   see subroutine mym_initialize
!       dfm(nx,nz,ny) : Diffusivity coefficient for momentum,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfh(nx,nz,ny) : Diffusivity coefficient for heat,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfq(nx,nz,ny) : Diffusivity coefficient for q^2,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       tcd(nx,nz,ny)   : Countergradient diffusion term for Theta_l
!                                                                  (K/s)
!       qcd(nx,nz,ny)   : Countergradient diffusion term for Q_w
!                                                              (kg/kg s)
!       pd?(nx,nz,ny) : Half of the production terms
!
!       Only tcd and qcd are defined at the center of the grid boxes
!
!     # DO NOT forget that tcd and qcd are added on the right-hand side
!       of the equations for Theta_l and Q_w, respectively.
!
!     Work arrays:        see subroutine mym_initialize and level2
!
!     # dtl, dqw, dtv, gm and gh are allowed to share storage units with
!       dfm, dfh, dfq, tcd and qcd, respectively, for saving memory.
!
*/

void MODULE_BL_MYNN_mp_MYM_TURBULENCE(
    int *, int *, int *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, int *, float *, float *, float *, int *, float *,
    float *, int *);

/*
! ==================================================================
!     SUBROUTINE  mym_predict:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       qke(nx,nz,ny) : qke at (n)th time level
!       tsq, ...cov     : ditto
!
!     Output variables:
!       qke(nx,nz,ny) : qke at (n+1)th time level
!       tsq, ...cov     : ditto
!
!     Work arrays:
!       qkw(nx,nz,ny)   : q at the center of the grid boxes        (m/s)
!       bp (nx,nz,ny)   : = 1/2*F,     see below
!       rp (nx,nz,ny)   : = P-1/2*F*Q, see below
!
!     # The equation for a turbulent quantity Q can be expressed as
!          dQ/dt + Ah + Av = Dh + Dv + P - F*Q,                      (1)
!       where A is the advection, D the diffusion, P the production,
!       F*Q the dissipation and h and v denote horizontal and vertical,
!       respectively. If Q is q^2, F is 2q/B_1L.
!       Using the Crank-Nicholson scheme for Av, Dv and F*Q, a finite
!       difference equation is written as
!          Q{n+1} - Q{n} = dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   )
!                        + dt/2*( Dv{n+1} - Av{n+1} - F*Q{n+1} ),    (2)
!       where n denotes the time level.
!       When the advection and diffusion terms are discretized as
!          dt/2*( Dv - Av ) = a(k)Q(k+1) - b(k)Q(k) + c(k)Q(k-1),    (3)
!       Eq.(2) can be rewritten as
!          - a(k)Q(k+1) + [ 1 + b(k) + dt/2*F ]Q(k) - c(k)Q(k-1)
!                 = Q{n} + dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   ),    (4)
!       where Q on the left-hand side is at (n+1)th time level.
!
!       In this subroutine, a(k), b(k) and c(k) are obtained from
!       subprogram coefvu and are passed to subprogram tinteg via
!       common. 1/2*F and P-1/2*F*Q are stored in bp and rp,
!       respectively. Subprogram tinteg solves Eq.(4).
!
!       Modify this subroutine according to your numerical integration
!       scheme (program).
!
!-------------------------------------------------------------------
*/

void MODULE_BL_MYNN_mp_MYM_PREDICT(int *, int *, int *, float *, float *,
                                   float *, float *, float *, float *, float *,
                                   float *, float *, float *, float *, float *,
                                   float *, float *, float *, float *, float *,
                                   float *, float *, int *);

/*
! ==================================================================
!     SUBROUTINE  mym_condensation:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       exner(nz)    : Perturbation of the Exner function    (J/kg K)
!                         defined on the walls of the grid boxes
!                         This is usually computed by integrating
!                         d(pi)/dz = h*g*tv/tref**2
!                         from the upper boundary, where tv is the
!                         virtual potential temperature minus tref.
!
!     Output variables:   see subroutine mym_initialize
!       cld(nx,nz,ny)   : Cloud fraction
!
!     Work arrays:
!       qmq(nx,nz,ny)   : Q_w-Q_{sl}, where Q_{sl} is the saturation
!                         specific humidity at T=Tl
!       alp(nx,nz,ny)   : Functions in the condensation process
!       bet(nx,nz,ny)   : ditto
!       sgm(nx,nz,ny)   : Combined standard deviation sigma_s
!                         multiplied by 2/alp
!
!     # qmq, alp, bet and sgm are allowed to share storage units with
!       any four of other work arrays for saving memory.
!
!     # Results are sensitive particularly to values of cp and rd.
!       Set these values to those adopted by you.
!
!-------------------------------------------------------------------
*/

void MODULE_BL_MYNN_mp_MYM_CONDENSATION(int *, int *, float *, float *, float *,
                                        float *, float *, float *, float *,
                                        float *, float *, float *, float *,
                                        int *, float *, float *, float *,
                                        float *, float *, float *, float *);

/*
       Subroutine MYNN_TENDENCIES.
*/

void MODULE_BL_MYNN_mp_MYNN_TENDENCIES(
    int *, int *, int *, int *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, int *, int *, int *, int *,
    float *, int *, int *, int *, int *);

/*
       SUBROUTINE mynn_mix_chem.
*/

void MODULE_BL_MYNN_mp_MYNN_MIX_CHEM(
    int *, int *, int *, int *, float *, float *, int *, int *, int *, int *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, int *);

/*
       Subroutine retrieve_exchange_coeffs.
*/

void MODULE_BL_MYNN_mp_RETRIEVE_EXCHANGE_COEFFS(int *, int *, float *, float *,
                                                float *, float *, float *,
                                                float *, float *);

/*
! ==================================================================
SUBROUTINE tridiag(n,a,b,c,d)

!! to solve system of linear eqs on tridiagonal matrix n times n
!! after Peaceman and Rachford, 1955
!! a,b,c,d - are vectors of order n
!! a,b,c - are coefficients on the LHS
!! d - is initially RHS on the output becomes a solution vector

!-------------------------------------------------------------------
*/

void MODULE_BL_MYNN_mp_TRIDIAG(int *, float *, float *, float *, float *);

/*
       SUBROUTINE mynn_bl_driver.
*/

void MODULE_BL_MYNN_mp_MYNN_BL_DRIVER(
    int *, int *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, int *, float *, float *, int *, int *, int *,
    int *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, int *, int *,
    float *, int *, int *, float *, float *, int *, int *, int *, int *, int *,
    int *, float *, float *, float *, float *, float *, float *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *);

/*
       SUBROUTINE mynn_bl_init_driver.
*/

void MODULE_BL_MYNN_mp_MYNN_BL_INIT_DRIVER(float *, float *, float *, float *,
                                           float *, float *, float *, float *,
                                           float *, int *, int *, int *, int *,
                                           int *, int *, int *, int *, int *,
                                           int *, int *, int *, int *, int *,
                                           int *, int *, int *, int *, int *,
                                           int *, int *);

/*
       SUBROUTINE get_pblh.
       !---------------------------------------------------------------
       !             NOTES ON THE PBLH FORMULATION
       !
       !The 1.5-theta-increase method defines PBL heights as the level at
       !which the potential temperature first exceeds the minimum potential
       !temperature within the boundary layer by 1.5 K. When applied to
       !observed temperatures, this method has been shown to produce PBL-
       !height estimates that are unbiased relative to profiler-based
       !estimates (Nielsen-Gammon et al. 2008). However, their study did not
       !include LLJs. Banta and Pichugina (2008) show that a TKE-based
       !threshold is a good estimate of the PBL height in LLJs. Therefore,
       !a hybrid definition is implemented that uses both methods, weighting
       !the TKE-method more during stable conditions (PBLH < 400 m).
       !A variable tke threshold (TKEeps) is used since no hard-wired
       !value could be found to work best in all conditions.
       !---------------------------------------------------------------
*/

void MODULE_BL_MYNN_mp_GET_PBLH(int *, int *, float *, float *, float *,
                                float *, float *, float *, float *);

/*
       SUBROUTINE StEM_mf
       ! ==================================================================
       ! Much thanks to Kay Suslj of NASA-JPL for contributing the original
   version ! of this mass-flux scheme. Considerable changes have been made from
   it's ! original form. Some additions include: !  1) scale-aware tapering as
   dx -> 0 !  2) transport of TKE (extra namelist option) !  3)
   Chaboureau-Bechtold cloud fraction & coupling to radiation (when icloud_bl >
   0) !  4) some extra limits for numerical stability ! This scheme remains
   under development, so consider it experimental code.
       !
*/

void MODULE_BL_MYNN_mp_STEM_MF(
    int *, int *, float *, float *, float *, int *, int *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, int *,
    float *, float *, float *, float *, int *, int *, float *);

/*
       SUBROUTINE POISSON
*/

void MODULE_BL_MYNN_mp_POISSON(int *, int *, int *, int *, float *, float *);

/*
       SUBROUTINE INIT_RANDOM_SEED.
*/

void MODULE_BL_MYNN_mp_INIT_RANDOM_SEED();

/*
       SUBROUTINE RANDOM_POISSON
       !**********************************************************************
       !     Translated to Fortran 90 by Alan Miller from:  RANLIB
       !
       !     Library of Fortran Routines for Random Number Generation
       !
       !                    Compiled and Written by:
       !
       !                         Barry W. Brown
       !                          James Lovato
       !
       !             Department of Biomathematics, Box 237
       !             The University of Texas, M.D. Anderson Cancer Center
       !             1515 Holcombe Boulevard
       !             Houston, TX      77030
       !
       ! Generates a single random deviate from a Poisson distribution with mean
   mu.
*/

void MODULE_BL_MYNN_mp_RANDOM_POISSON(float *, int *, int *);

/*
       SUBROUTINE CONDENSATION_EDMF.
*/

void MODULE_BL_MYNN_mp_CONDENSATION_EDMF(float *, float *, float *, float *,
                                         float *, float *);

/*
       SUBROUTINE SCALE_AWARE.
*/

void MODULE_BL_MYNN_mp_SCALE_AWARE(float *, float *, float *, float *);

/*
! ===================================================================
! ===================================================================
! This is the mass flux part of the TEMF scheme from module_bl_temf.F,
! adapted for the MYNN context by Wayne Angevine June 2015.
! Variable strategy:  TEMF external variables that have semantically
! comfortable counterparts in the MYNN-EDMF context have been changed to
! use those names.  Otherwise the TEMF variable names have been kept but
! redefined as local variables.  Only "moist" vars are used, whether
! updraft condenses or not.  Some former local vars are replaced with
! externals.
!
! (Partial) list of conversions:
!    wupd_temfx -> moist_w
!    thup_temfx -> moist_thl
!    qtup_temfx -> moist_qt
!    qlup_temfx -> moist_qc
!    cf3d_temfx -> cldfra_bl1d
!    au -> moist_a
*/

void MODULE_BL_MYNN_mp_TEMF_MF(
    int *, int *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, float *, float *, float *, float *, float *,
    float *, float *, float *, int *, float *, float *, float *, float *, int *,
    int *, float *, int *, int *, int *, int *, int *, int *);

/*
 */
}

#endif

#endif /*__GMS_MODULE_BL_MYNN_C_IFACE_H__*/