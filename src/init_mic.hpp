/***
  *  Project:
  *
  *  File: init_mic.cuh
  *  Created: Apt 02, 2013
  *  Modified: Sat 06 Apr 2013 11:24:47 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __INIT_MIC_HPP__
#define __INIT_MIC_HPP__

namespace hig {

	void init_mic() {
		#pragma offload_transfer target(mic:0)
	} // init_mic()

} // namespace hig

#endif // __INIT_MIC_HPP__
