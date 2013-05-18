/***
  *  Project:
  *
  *  File: fit_param_map.hpp
  *  Created: May 17, 2013
  *  Modified: Fri 17 May 2013 09:45:44 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __FIT_PARAM_MAP__
#define __FIT_PARAM_MAP__

#include <unordered_map>
#include <string>

#include "hig_input.hpp"


namespace hig {

	class FitParamMap {		// make it singleton? ...
		private:
			std::unordered_map <std::string, float_t> mapper_;
			// TODO: handle different kinds of values ...
			// have different "mapper_"s for each type of value ... ?

		public:
			FitParamMap() { }
			~FitParamMap() { }

			// this map is created/populated during the input construction

			// accessing and assigning
			float_t& operator[](std::string str) { return mapper_[str]; }

			friend class HiGInput;
	}; // class FitParamMap

} // namespace hig

#endif // __FIT_PARAM_MAP__
