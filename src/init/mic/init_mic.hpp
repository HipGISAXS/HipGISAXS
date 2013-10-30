/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: init_mic.cuh
 *  Created: Apt 02, 2013
 *  Modified: Tue 16 Jul 2013 11:51:04 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __INIT_MIC_HPP__
#define __INIT_MIC_HPP__

namespace hig {

	void init_mic() {
		#pragma offload_transfer target(mic:0)
	} // init_mic()

} // namespace hig

#endif // __INIT_MIC_HPP__
