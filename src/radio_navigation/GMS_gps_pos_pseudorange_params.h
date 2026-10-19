

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
                                                                                             
};

}

}

#endif /*__GMS_GPS_POS_PSEUDORANGE_PARAMS_H__*/