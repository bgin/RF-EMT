#ifndef _ATM_ENUMERATION_H
#define _ATM_ENUMERATION_H
/*******************************************************************************
 * ALMA - Atacama Large Millimiter Array
 * (c) Instituto de Estructura de la Materia, 2009
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * "@(#) $Id: ATMEnumerations.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error "This is a C++ include file and cannot be used from plain C"
#endif

#include "ATMCommon.h"

ATM_NAMESPACE_BEGIN

/** \file
 *  \brief Global enumeration objects
 *
 *  The ASDM uses the the item netSideband to describe the sidebands. It appears
 *  more convenient to use a pair of items, the side of the sideband and its type.
 *  The following two enumerations describes this.
 *
 */

/** \brief Defines the side of a sideband.
 *
 *  Heterodyne receivers have two sidebands, a lower and an upper sideband. On the
 *  other hand the detectors do not have sidebands. The Science Data Model assign
 *  one spectral window per sideband. Hence, with detectors, there is a single
 *  spectral window. SidebandSide differentiates between these various cases.<br>
 *
 *  The codes have the following meaning:
 *    <ul>
 *     <li> 0  No sideband (eg for a detector)
 *     <li> 1  Lower sideband (LSB)
 *     <li> 2  Upper sideband (USB)
 *    </ul>
 */
enum SidebandSide{ NOSB=0,         //!< No sideband (eg for HEMT)
                   LSB=1,          //!< Lower sideband (LSB)
                   USB=2           //!< Upper sideband (USB)
};


/** \brief Defines the kind of sideband separations.
 *
 *  With heterodyne receivers there are always two sidebands. In the
 *  interferometric case it is possible to separate these two sidebands.
 *  When there is no sideband separation two types are considered,
 *  single sideband (SSB) when the image sideband has been rejected
 *  and double sideband (DSB) when the siband gain ratio is about
 *  half half. In the case of DSB only one sideband is retained
 *  for the data because the other sideband has identical data.
 *  When modeling the atmospheric contribution both sidebands need to be
 *  considered because the contribution is not the same in the two
 *  sidebands. <br>
 *
 *  The codes have the following meaning:
 *   <ul>
 *     <li> -1 No type (NOTYPE): type is not relevant with a SidebandSide equal to NOSB
 *     <li> 0  Double sidebands (DSB): Case with no sideband rejection
 *     <li> 1  Single sideband (SSB): Case with sideband rejection
 *     <li> 2  Two sideband (TWOSB): Case with sideband separation (only for interferometric modes)
 *    </ul>
 *  </ul>
 */
enum SidebandType{ NOTYPE=-1,      //!< Not relevant (no sideband)
                   DSB=0,          //!< Double sideband (no rejection)
                   SSB=1,          //!< Single sideband (rejection)
                   TWOSB=2         //!< two sidebands (sideband separation)
};

 
/** \brief Defines the type of atmosphere
 *
 *  
 *
 *  The codes have the following meaning:
 *   <ul>
 *     <li> 1  Tropical (tropical)
 *     <li> 2  Mid-latitude Summer (midlatSummer) 
 *     <li> 3  Mid-latitude Winter (midlatWinter) 
 *     <li> 4  Sub-arctic Summer (subarcticSummer) 
 *     <li> 5  Sub-arctic Winter (subarcticWinter) 
 *     <li> 6  typeATM_end
 *    </ul>
 *  </ul>
 */

enum typeAtm_t { tropical=1, 
 		 midlatSummer, 
 		 midlatWinter, 
 		 subarcticSummer, 
 		 subarcticWinter, 
 		 typeATM_end 
};


ATM_NAMESPACE_END

#endif /*!_ATM_ENUMERATION_H*/
