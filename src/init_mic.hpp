/***
  *  Project:
  *
  *  File: init_mic.cuh
  *  Created: Apt 02, 2013
  *  Modified: Tue 02 Apr 2013 03:37:52 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __INIT_MIC_HPP__
#define __INIT_MIC_HPP__

namespace hig {

	void init_mic() {
		std::cout << "-- Waking up MIC(s) ..." << std::flush << std::flush;
		#pragma offload_transfer target(mic:0)
		std::cout << " done." << std::endl;
	} // init_mic()

} // namespace hig

#endif // __INIT_MIC_HPP__
