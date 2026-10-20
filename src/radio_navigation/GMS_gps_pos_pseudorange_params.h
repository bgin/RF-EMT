

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

#ifndef __GMS_GPS_POS_PSEUDORANGE_PARAMS_H__
#define __GMS_GPS_POS_PSEUDORANGE_PARAMS_H__ 191020261346

#include <string>
#include <tuple>
#include <cstdint>

namespace file_info 
{

     static const unsigned int GMS_GPS_POS_PSEUDORANGE_PARAMS_MAJOR = 1;
     static const unsigned int GMS_GPS_POS_PSEUDORANGE_PARAMS_MINOR = 1;
     static const unsigned int GMS_GPS_POS_PSEUDORANGE_PARAMS_MICRO = 0;
     static const unsigned int GMS_GPS_POS_PSEUDORANGE_PARAMS_FULLVER =
       1000U*GMS_GPS_POS_PSEUDORANGE_PARAMS_MAJOR+100U*GMS_GPS_POS_PSEUDORANGE_PARAMS_MINOR+
       10U*GMS_GPS_POS_PSEUDORANGE_PARAMS_MICRO;
     static const char GMS_GPS_POS_PSEUDORANGE_PARAMS_CREATION_DATE[] = "19-10-2026 01:52PM +00200 (MON 19 OCT 2026 GMT+2)";
     static const char GMS_GPS_POS_PSEUDORANGE_PARAMS_BUILD_DATE[]    = __DATE__; 
     static const char GMS_GPS_POS_PSEUDORANGE_PARAMS_BUILD_TIME[]    = __TIME__;
     static const char GMS_GPS_POS_PSEUDORANGE_PARAMS_AUTHOR[]        = "Programmer: Bernard Gingold, beniekg@gmail.com";
     static const char GMS_GPS_POS_PSEUDORANGE_PARAMS_SYNOPSIS[]      = "Based upon the chapter 2 (part 1) of the book titled: Global Positioning System,Volume 1, Theory and Applications, \
                                                                         edited by: James J. Spilker, Bradford W. Parkinson";
}

namespace gms 
{

namespace radio_navigation
{

alignas(64) struct GPS_pos_pseudorange_params_t 
{
    using tuple_1ui_2s = std::tuple<std::uint32_t,std::string,std::string>;

    const tuple_1ui_2s m_sat_pos      = std::make_tuple(1U,"xsi or rsi", "Position of satellite i in either x, y, z or radius vector coordinates,\
                                                                       respectively. The position of the satellites at the time of \
                                                                       transmission is xsi (tsi) or rsi (tsi).");
    const tuple_1ui_2s m_user_pos     = std::make_tuple(2U,"xu or ru",   "User position in either x, y, z or radius vector coordinates. User \
                                                                       position may also vary with time. User position at time of \
                                                                       reception is xui (tRi) or rui (tRi).");
    const tuple_1ui_2s m_sat_tx_time  = std::make_tuple(3U,"tsi","True time of transmission from satellite i. This parameter is the \
                                                                  true time of transmission and may not be exactly the same \
                                                                  as indicated by the satellite.");
    const tuple_1ui_2s m_sat_clk_read = std::make_tuple(4U,"t'si","Actual satellite clock reading");
    const tuple_1ui_2s m_sat_clk_bias = std::make_tuple(5U,"delta bi","Satellite bias clock error (expressed in s)");      
    const tuple_1ui_2s m_gps_sys_time = std::make_tuple(6U,"t","GPS system time");
    const tuple_1ui_2s m_act_user_time= std::make_tuple(7U,"t'si", "Actual user clock time at time of reception of signal");
    const tuple_1ui_2s m_true_user_time=std::make_yuple(8U,"tu", "True user time at time of reception");
    const tuple_1ui_2s m_user_clk_bias= std::make_tuple(9U,"bu","User clock bias—can vary with time");
    const tuple_1ui_2s m_true_val_psr = std::make_tuple(10U,"rho_iT", "True value of pseudorange");
    const tuple_1ui_2s m_meas_err_psr = std::make_tuple(11U,"rho_i", "Measured pseudorange with various error contributors");
    const tuple_1ui_2s m_rx_thermal_n = std::make_tuple(12U,"n(t)", "Receiver thermal noise");
    const tuple_1ui_2s m_sat_bias_clk_err = std::make_tuple(13U,"delta_bi", "Satellite bias clock error (expressed in s)");
    const tuple_1ui_2s m_sat_bias_pos_err = std::make_tuple(14U,"delta_Di", "Satellite position bias error effect on range");
    const tuple_1ui_2s m_rx_psr_noise_err = std::make_tuple(15U,"mu_i", "Receiver pseudorange measurement noise error for satellite i,s");
    const tuple_1ui_2s m_ionosph_delay    = std::make_tuple(16U,"delta_Ii", "Ionospheric delay (expressed in s)");
    const tuple_1ui_2s m_troposph_delay   = std::make_tuple(17U,"delta_Ti", "Tropospheric delay (expressed in s)");
    const tuple_1ui_2s m_rel_time_correct = std::make_tuple(18U,"delta_mui","Relativistic time correction (expressed in s)");
    const tuple_1ui_2s m_velocity_light   = std::make_tuple(18U,"c","Velocity of Light");
    const tuple_1ui_2s m_DOP              = std::make_tuple(19U,"DOP", "Dilution of Precession");
    const tuple_1ui_2s m_VDOP             = std::make_tuple(20U,"VDOP", "Vertical geometric dilutions of precession");
    const tuple_1ui_2s m_HDOP             = std::make_tuple(21U,"HDOP", "Horizontal geometric dilutions of precession"); 
    const tuple_1ui_2s m_PDOP             = std::make_tuple(22U,"PDOP", "Position geometric dilutions of precession");  
    const tuple_1ui_2s m_GDOP             = std::make_tuple(23U,"GDOP", "Geometric dilution of precession includes both position and time error effects");
    const tuple_1ui_2s m_meas_carr_ph_off = std::make_tuple(24U,"phi_i", "Measured carrier phase offset as received from satellite i");
    const tuple_1ui_2s m_true_rx_cph_off  = std::make_tuple(25U,"phi_iT","True received carrier phase offset, at frequency Omega0");
    const tuple_1ui_2s m_gps_carr_wavel   = std::make_tuple(26U,"gamma", "Wavelength of GPS carriers");      
    const tuple_1ui_2s m_carr_cyc_ambig   = std::make_tuple(27U,"n", "Carrier phase count cycles ambiguity");
    const tuple_1ui_2s m_true_user_rng    = std::make_tuple(28U,"D_i", "True range from satellite i to user");  
    inline constexpr GPS_pos_pseudorange_params_t & operator=(const GPS_pos_pseudorange_params_t &) /* noexcept */ = delete;
    inline constexpr GPS_pos_pseudorange_params_t & operator=(GPS_pos_pseudorange_params_t &&) /* noexcept */ = delete;                                                                        
};

}

}

#endif /*__GMS_GPS_POS_PSEUDORANGE_PARAMS_H__*/