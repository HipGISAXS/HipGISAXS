/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: enums.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:01 AM PST
 *  Description: Defines the common enums used.
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _ENUMS_HPP_
#define _ENUMS_HPP_

namespace hig{

  enum AlgoType{
    fit_brute_force=0,
    fit_cg,
    fit_pounders,
    fit_genetic,
    fit_swarm,
    stat_baysian,
    stat_rmc,
    ml_cluster,
    ml_classify
  };

  enum ImgMode{
    q_space=0,
    angles,
    pixels
  };


} /* namespace hig */

#endif /* ENUMS_HPP_ */
