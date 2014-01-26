/**
 *  Project: WOO Timer Library
 *
 *  File: wootimers.hpp
 *  Created: Nov 21, 2012
 *  Modified: Wed 17 Jul 2013 10:26:44 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Copyright (c) 2012-2013 Abhinav Sarje
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE file.
 */

#ifndef _WOOTIMERS_HPP_
#define _WOOTIMERS_HPP_

namespace woo {

class WooTimer {		// an abstract class
	protected:
		double start_;
		double stop_;
		double elapsed_;	// in lowest resolution of respective timers used
		bool is_running_;

	public:
		virtual ~WooTimer() { }

		virtual void reset() = 0;

		virtual void start() = 0;
		virtual void stop() = 0;
		virtual double lap() = 0;			// lap in lowest resolution

		virtual double elapsed_sec() = 0;	// in seconds		10^0
		virtual double elapsed_msec() = 0;	// in miliseconds	10^3
		virtual double elapsed_usec() = 0;	// in microseconds	10^6
		virtual double elapsed_nsec() = 0;	// in nanoseconds	10^9
}; // class WooTimer

} // namespace woo


#endif // _WOOTIMERS_HPP_
