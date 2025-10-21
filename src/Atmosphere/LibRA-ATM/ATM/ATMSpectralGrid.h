#ifndef _ATM_SPECTRALGRID_H
#define _ATM_SPECTRALGRID_H
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
 * "@(#) $Id: ATMSpectralGrid.h Exp $"
 *
 * who       when      what
 * --------  --------  ----------------------------------------------
 * pardo     24/03/09  created
 */

#ifndef __cplusplus
#error "This is a C++ include file and cannot be used from plain C"
#endif

#include "ATMCommon.h"
#include "ATMEnumerations.h"
#include "ATMFrequency.h"
#include <string>
#include <vector>

ATM_NAMESPACE_BEGIN
/*!  \brief Spectral grid defined as a list of frequencies with basic complementary
 *          informations.
 *
 *  This grid is composed by a set frequency channels (or pixels) which may or may not
 *  be at a regular interval. The grid is composed of one or more spectral windows.
 *  A spectral window is defined as a group of channels, this group having the following
 *  attributes:
 *  <ul>
 *  <li> a number of channels
 *  <li> a reference channel
 *  <li> the frequency for that reference channel
 *  <li> a channel separation (it is set to 0 if there is only one channel or when the
 *      channels are at irregular intervals.
 *  <li> the intermediate frequency (if NetSideband is greater than -1)
 *  To every spectral window is assigned an identifier (zero-based) which can be used
 *  when retrieving the parameters of a given window.<br>
 *  In case of spectral windows with regularily spaced channels this reference frequency
 *  can be associated to a position in channel number units.
 *
 *  There are several different constructors which can be used to define a spectral
 *  window. One of these must be used to construct the first spectral window. Would
 *  more than a single spectral be necessary, they are appended with the add methods.
 *  When there are two sidebands, every of these corresponds to a single spectral window.
 *  Note that these two sidebands can be constructed using a single method if the intermediate
 *  frequency is an input parameter.
 *  \note The Alma Science Data Model (ASDM) spectral windows derive from spectral grids,
 *  with extra informations such as the channel width, the resolution the noise bandwidth
 *  etc...
 *
 */
class SpectralGrid
{
public:
  SpectralGrid(const Frequency &oneFreq);

  /** A full constructor to be used in cases the channel separation is uniform. The
   *  position of the reference frequency is set at the center of the specified
   *  reference channel. This reference channel does not need to be necessarily
   *  in the range [1,numChan].
   * @pre no spectral window has yet been defined
   * @param numChan Number of channels for the spectral window
   * @param refChan Reference channel for the spectral window
   * @param refFreq Frequency for this reference channel
   * @param chanSep Frequency increment between two adjacent channels
   * @post  the spectral window has been defined; it is taken as a spectral window
   *        with no sideband. Its spectral window identifier spwId is 0.
   */
  SpectralGrid(unsigned int numChan,
               unsigned int refChan,
               const Frequency &refFreq,
               const Frequency &chanSep);

  /** A full constructor to be used in cases the channel separation is not uniform.
   *  In this case the reference frequency may not coincide with the center
   *  of a channel. The frequency units (freqUnits) is common for the reference
   *  frequency (refFReq) and for the numChan frequencies in the list chanFreq.
   *  Known units are GHz, MHz, kHz and Hz.
   * @pre no spectral window has yet been defined
   * @param numChan Number of channels for the spectral window
   * @param refFreq The frequency of the reference channel
   * @param chanFreq The frequencies of the individual channels
   * @param freqUnits The frequency units used to define refFreq and chanFreq
   * @post  the spectral window has been defined; it is taken as a spectral window
   *        with no sideband. Its spectral window identifier spwId is 0.
   */
  SpectralGrid(unsigned int numChan,
               double refFreq,
               double* chanFreq,
               Frequency::Units freqUnits);
  SpectralGrid(double refFreq, const std::vector<double> &chanFreq, Frequency::Units freqUnits);
  SpectralGrid(const std::vector<double> &chanFreq, Frequency::Units freqUnits);
  SpectralGrid(const std::vector<Frequency> &chanFreq);

  /** A full constructor to be used in cases of two sibands (separated or not)
   *  The reference frequency, may not coincide with the center
   *  of a channel, is the frequency used to track on.
   *  The frequency units (freqUnits) is common for the reference
   *  frequency (refFReq) and for the numChan frequencies in the list chanFreq.
   *  Known units are GHz, MHz, kHz and Hz.
   * @pre no spectral window has yet been defined
   * @param numChan Number of channels for the spectral window
   * @param refFreq Frequency of the reference channel
   * @param chanSep Frequency interval between two adajacent channels
   * @param intermediateFreq Intermediate Frequency
   * @param SidebandSide code (see #SidebandSide for details)
   * @param SidebandType code (see #SidebandType for details)
   * @post  <b>Two</b> spectral windows have been defined, one for each sideband
   *        Their spectral window identifiers are 0 and 1. The identifier is 0 and 1 for the lower
   *        and upper sidebands if the input netSideband is 0, 1, 3 or 5 else it is 1 and 0 respectively.
   */
  SpectralGrid(unsigned int numChan,
               unsigned int refChan,
               const Frequency &refFreq,
               const Frequency &chanSep,
               const Frequency &intermediateFreq,
               const SidebandSide &sbSide,
               const SidebandType &sbType);

  /** A full constructor to be used in cases the channel separation is not uniform.
   *  In this case the reference frequency coincides with the center of a channel,
   *  the one at the refChan position in the numChan elements of the list of
   *  frequencies chanFreq. Hence refChan must be in the range [1,numChan]. The
   *  frequency units (freqUnits) is common for the reference
   *  frequency (refFReq) and for the numChan frequencies in the list chanFreq.
   *  Known units are GHz, MHz, kHz and Hz.
   * @pre no spectral window has yet been defined
   * @param numChan Number of channels for the spectral window
   * @param refChan Reference channel for the spectral window
   * @param chanFreq The frequencies of the individual channels
   * @param freqUnits The frequency units used to define refFreq and chanFreq
   * @post  the spectral window has been defined; it is taken as a spectral window
   *        with no sideband. Its spectral window identifier spwId is 0.
   */
  SpectralGrid(unsigned int numChan,
               unsigned int refChan,
               double* chanFreq,
               Frequency::Units freqUnits);

  SpectralGrid();

  SpectralGrid(const SpectralGrid &);

  ~SpectralGrid();

  /** Add a new spectral window, uniformly sampled, this spectral window having no sideband.
   * @pre   at least one spectral window has already been defined
   * @param numChan Number of channels for the spectral window
   * @param refChan Reference channel for the spectral window
   * @param refFreq Frequency for this reference channel
   * @param chanSep Frequency increment between two adjacent channels
   * @return the identifier for this spectral window
   * @post  the spectral window has been defined and appended to the set of spectral windows;
   *        this new spectral window is taken as one with no sideband.
   */
  unsigned int add(unsigned int numChan,
                   unsigned int refChan,
                   const Frequency &refFreq,
                   const Frequency &chanSep);

  /** Add a new spectral window, this spectral window having no sideband.
   * @pre   at least one spectral window has already been defined
   * @param numChan Number of channels for the spectral window
   * @param refChan Reference channel for the spectral window
   * @param chanFreq The frequency for every channel
   * @param freqUnits The frequency units (GHz, MHz or kHz, else assumed in Hz) used for chanFreq
   * @return the identifier for this spectral window
   * @post  the spectral window has been defined and appended to the set of spectral windows;
   *        this new spectral window is taken as one with no sideband.
   */
  unsigned int add(unsigned int numChan,
                   unsigned int refChan,
                   double* chanFreq,
                   Frequency::Units freqUnits);

  /** Add a new spectral window, this spectral window having no sideband.
   * @pre   at least one spectral window has already been defined
   * @param numChan Number of channels for the spectral window
   * @param refFreq The frequency at the reference channel for the spectral window
   * @param chanFreq The frequency for every channel
   * @param freqUnits The frequency units (GHz, MHz or kHz, else assumed in Hz) used for refFreq and chanFreq
   * @return the identifier for this spectral window
   * @post  the spectral window has been defined and appended to the set of spectral windows;
   *        this new spectral window is taken as one with no sideband.
   */
  unsigned int add(unsigned int numChan,
                   double refFreq,
                   double* chanFreq,
                   Frequency::Units freqUnits);
  unsigned int add(unsigned int numChan,
                   double refFreq,
                   const std::vector<double> &chanFreq,
                   Frequency::Units freqUnits);
  unsigned int add(const std::vector<Frequency> &chanFreq)
  {
    unsigned int spwid;
    std::vector<double> v;
    for(unsigned int i = 0; i < chanFreq.size(); i++) {
      v.push_back(chanFreq[i].get(Frequency::UnitGigaHertz));
    }
    spwid = add(chanFreq.size(), chanFreq[0].get(Frequency::UnitGigaHertz), v, Frequency::UnitGigaHertz);
    return spwid;
  }
  /** Add two new spectral windows, one spectral window per sideband.
   * @pre   at least one spectral window has already been defined
   * @param numChan Number of channels for the spectral window in one sideband
   * @param refFreq The frequency at the reference channel for that spectral window
   * @param chanSep The channel separation (unsigned)
   * @param intermediateFreq the Intermediate Frequency
   * @return the netSideband for this first spectral window (the second spectral will have the corresponding
   *         netSideband for the other sideband.
   * @post  the spectral window has been defined and appended to the set of spectral windows;
   *        this new spectral window is taken as one with no sideband.
   * \note For SSB receivers the first sideband is the signal sideband and in this case refFreq is usually
   *       the frequency for the middle of that sideband. For DSB refFreq may be in between the two sideband.

   */
  void add(unsigned int numChan,
           unsigned int refChan,
           const Frequency &refFreq,
           const Frequency &chanSep,
           const Frequency &intermediateFreq,
           const SidebandSide &sbSide,
           const SidebandType &sbType);

  //@{
  /** Accessor to the number of spectral wondows
   * @return number of spectral windows
   */
  unsigned int getNumSpectralWindow() const;

  /** Accessor to the number of frequency points for the first spectral window
   * @return number of frequency channels
   */
  unsigned int getNumChan() const;
  /** Accessor to the number of frequency points for a given spectral window
   * @param spwId spectral window identifier (0-based)
   * @return number of frequency channels
   */
  unsigned int getNumChan(unsigned int spwId) const;

  /** Accessor to the reference channel of the first spectral window
   * @return the reference channel
   * \note with the constructor SpectralGrid( int numChan, double refFreq, double* chanFreq, std::string freqUnits)
   *       there is no way to determine the reference channel if the grid is not regularily sampled! would that be
   *       the case, the returned value is 0.
   */
  unsigned int getRefChan() const;
  /** Accessor to the reference channel for a given spectral window
   * @param spwId spectral window identifier (0-based)
   * @return the reference channel
   */
  unsigned int getRefChan(unsigned int spwId) const;

  // Frequency getRefFreq();
  /** Accessor to the reference frequency (Hz) for the first spectral window
   * @return the frequency at the reference channel position
   */
  Frequency getRefFreq() const;

  /** Accessor to the reference frequency (Hz) for a given spectral window
   * @param spwId spectral window identifier (0-based)
   * @return the frequency at the reference channel position
   */
  Frequency getRefFreq(unsigned int spwId) const;

  /** Accessor to the channel separation for regularily spaced grids (for the first spectral window)
   * @return the channel separation (Hz)
   */
  Frequency getChanSep() const;

  /** Accessor to the channel separation for regularily spaced grids (for a given spectral window)
   * @param spwId spectral window identifier (0-based)
   * @return the channel separation (Hz)
   */
  Frequency getChanSep(unsigned int spwId) const;

  /** Accessor to the frequency (Hz) for a given grid point (for the first spectral window)
   * @param chanNum    the channel number (grid units)
   * @ return the frequency (Hz) corresponding to the center of the channel
   */
  Frequency getChanFreq(unsigned int chanNum) const;
  Frequency getChanWidth(unsigned int chanNum) const;

  /** Accessor to the frequency (Hz) for a given grid point for the specified spectral window
   * @param spwId     spectral window identifier (0-based)
   * @param chanNum   the channel number (grid units)
   * @ return the frequency (Hz) corresponding to the center of the channel
   */
  Frequency getChanFreq(unsigned int spwId, unsigned int chanNum) const;
  Frequency getChanWidth(unsigned int spwId, unsigned int chanNum) const;


  /** Accessor to the frequencies in the specified units for a given channel index (0-based) for the
   * the specified spectral window its corresponding other sideband.
   * @param spwId     spectral window identifier (0-based)
   * @param chanNum   the channel number (grid units)
   * @param freqUnits the requested units
   * @return  the frequencies corresponding to the center of the channel in the specified units for
   *          the input spectral window and for the corresponding channel in the other sideband
   */
  std::vector<double> getSbChanFreq(unsigned int spwId,
                               unsigned int chanNum,
                               const std::string &freqUnits) const;

  /** Accessor to retrieve the spectral grid of a spectral window
   * @param spwId     spectral window identifier (0-based)
   * @return a std::vector of numChan frequencies (Hz)
   */
  std::vector<double> getSpectralWindow(unsigned int spwId) const;

  /** Method to get the grid position for a given frequency specified in Hz (the first spectral window)
   * @return the grid position
   */
  double getChanNum(double freq) const;
  /** Method to get the grid position for a given frequency specified in Hz for the specified spectral window
   * @param freq the frequency (Hz) for the grid position
   * @param spwId spectral window identifier (0-based)
   * @return the grid position
   */
  double getChanNum(unsigned int spwId, double freq) const;

  /** Method to get the frequency range encompassing the list of frequency grid points (for the first spectral window)
   * \note In case of irregular sampling the return value is the difference between the frequencies of the channels
   *       with the highest one and the channel with the lowedt one else it the product of the number of channels
   *       times the channel separation
   * @return the frequency bandwidth (Hz)
   */
  Frequency getBandwidth() const;

  /** Method to get the frequency range encompassing the list of frequency grid points for the specified spectral window
   * \note In case of irregular sampling the return value is the difference between the frequencies of the channels
   *       with the highest one and the channel with the lowedt one else it the product of the number of channels
   *       times the channel separation
   * @param spwId spectral window identifier (0-based)
   * @return the frequency bandwidth (Hz)
   */
  Frequency getBandwidth(unsigned int spwId) const;

  /** Method to get the frequency (Hz) for the point at the lowest frequency (for the first spectral window)
   * @return the frequency of the channel at the lowest frequency
   */
  Frequency getMinFreq() const;

  /** Method to get the frequency in the specified units for the channel at the lowest frequency for the
   * specified spectral window
   * @param spwId spectral window identifier (0-based)
   * @return the frequency (Hz) of the channel at the lowest frequency
   */
  Frequency getMinFreq(unsigned int spwId) const;

  /** Method to get the frequency (Hz) for the point at the largest frequency (for the first spectral window)
   * @return the frequency (Hz) of the channel at the highest frequency
   */
  Frequency getMaxFreq() const;
  /** Method to get the frequency in the specified units for the point at the largest frequency (for the first
   * spectral window)
   * @param  units the requested units
   * @return the frequency (Hz) of the channel at the highest frequency
   */
  Frequency getMaxFreq(unsigned int spwId) const;

  /** Method to know if the spectral grid is regular or not (the first spectral window)
   * @return true if uniformly sampled, else false
   */
  bool isRegular() const;
  /** Method to know if the spectral grid is regular or not for the specified spectral window
   * @param spwId spectral window identifier (0-based)
   * @return true if uniformly sampled, else false
   */
  bool isRegular(unsigned int spwId) const;

  /** Accessor for the side of the sideband
   * @param spwId spectral window identifier (0-based)
   * @return the side of the sideband
   * \note Possible result is no sideband (NoSB) or lower sideband (LSB) or upper sideband (USB)
   */
  std::string getSidebandSide(unsigned int spwId) const;

  /** Accessor to the nature(s) of the associated spectral window(s)
   * @pre the spectral window must have an associated sideband. Would that not be the
   *      case the returned std::vector has a size of 0.
   * @param spwId spectral window identifier (0-based)
   * @return the associated nature(s) of the associated spectral windows
   */
  std::vector<std::string> getAssocNature(unsigned int spwId) const;

  /** Accessor to the identifier of the associated spectral window(s)
   * @pre the spectral window must have an associated spectral window. Would that not be the
   *      case the returned std::vector has a size of 0.
   * @param spwId spectral window identifier (0-based)
   * @return the identifiers of the its associated spectral windows
   */
  std::vector<unsigned int> getAssocSpwId(unsigned int spwId) const;

  std::vector<unsigned int> getAssocSpwIds(const std::vector<unsigned int> &spwIds) const;

  /** Accessor for the type of  sideband
   * @pre the spectral window must have an associated spectral window. Would that not be the
   *      case the returned std::vector has a size of 0.
   * @param spwId spectral window identifier (0-based)
   * @return the type of the sideband
   * \note Possible result is double sideband (DSB) or single sideband (SSB) or two sidebands (2SB).
   * 2SB implies sideband separation which is possible only in the interferometric case
   */
  std::string getSidebandType(unsigned int spwId) const;

  /** Accessor for the side of the sideband and its type
   * @pre the spectral window must have a sideband side and a sideband type. Would that not be the
   *      case the returned std::string has a size of 0.
   * @param spwId spectral window identifier (0-based)
   * @return the side and the type of the sideband
   */
  std::string getSideband(unsigned int spwId) const;

  double getLoFrequency() const;
  double getLoFrequency(unsigned int spwId) const;

  //@}

  bool operator ==(const SpectralGrid&) const;

protected:
  Frequency::Units freqUnits_; //!< The frequency inits (always Frequency::UnitHertz)
  std::vector<double> v_chanFreq_; //!< Channel frequencies of ALL the channels (i.e. all spectral window appended)

  std::vector<unsigned int> v_numChan_; //!< number of channels for every spectral window
  std::vector<unsigned int> v_refChan_; //!< reference channel for every spectral window
  std::vector<double> v_refFreq_; //!< frequency at reference channel for every spectral window
  std::vector<double> v_chanSep_; //!< channel separation for every spectral window
  std::vector<double> v_maxFreq_; //!< frequency maximum for every spectral window
  std::vector<double> v_minFreq_; //!< frequency minimum for every spectral window
  std::vector<double> v_intermediateFrequency_; //<! intermediate frequency of the band for every spectral window
  std::vector<double> v_loFreq_; //<! LO frequency

  std::vector<SidebandSide> v_sidebandSide_; //<!  NOSB=0, LSB=1, USB=2
  std::vector<SidebandType> v_sidebandType_; //<!  DSB=0, SSB=1, TWOSB=2
  std::vector<std::vector<unsigned int> > vv_assocSpwId_; //<!  associated spectral window Id(s)
  std::vector<std::vector<std::string> > vv_assocNature_; //<!  corresponding associated nature

  std::vector<unsigned int> v_transfertId_;

private:

  void appendChanFreq(unsigned int numChan, double* chanFreq);
  void appendChanFreq(unsigned int numChan, const std::vector<double> &chanFreq);
  bool wrongSpwId(unsigned int spwId) const;
};

inline bool SpectralGrid::operator ==(const SpectralGrid & a) const
{
  bool equals = false;

  if(freqUnits_ != a.freqUnits_) return equals;
  if(v_refChan_ != a.v_refChan_) return equals;
  if(v_chanSep_ != a.v_chanSep_) return equals;
  if(v_numChan_ != a.v_numChan_) return equals;
  if(v_maxFreq_ != a.v_maxFreq_) return equals;
  if(v_minFreq_ != a.v_minFreq_) return equals;
  if(!equals)

  if(v_chanFreq_ != a.v_chanFreq_) return equals;

  equals = true;
  return equals;
}; // class SpectralGrid

ATM_NAMESPACE_END

#endif /*!_ATM_SPECTRALGRID_H*/
