/***
  *  Project: WOO Timer Library
  *
  *  File: woo_boostchronotimers.hpp
  *  Created: Nov 21, 2012
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "wootimers.hpp"
#include <boost/chrono.hpp>

#ifndef _WOOBOOSTCHRONOTIMERS_HPP_
#define _WOOBOOSTCHRONOTIMERS_HPP_

namespace woo {

class BoostChronoTimer : public WooTimer {
	// use Boost's chrono
	private:
		boost::chrono::duration<long long, boost::nano> chstart_;
		boost::chrono::duration<long long, boost::nano> chstop_;

		boost::chrono::steady_clock::time_point startpoint_;
		boost::chrono::steady_clock::time_point stoppoint_;

	public:
		BoostChronoTimer() {
			start_ = 0.0; stop_ = 0.0; elapsed_ = 0.0; is_running_ = false; is_paused_ = false;
		} // BoostChronoTimer()

		~BoostChronoTimer() { }

		void reset() {
			start_ = 0.0; stop_ = 0.0; elapsed_ = 0.0; is_running_ = false; is_paused_ = false;
		} // reset()

		void start() {
			reset();
			startpoint_ = boost::chrono::steady_clock::now();
			chstart_ = startpoint_.time_since_epoch();
			start_ = chstart_.count();
			is_running_ = true;
		} // start()

		void stop() {
			if(!is_running_) {
				std::cerr << "error: timer is not running" << std::endl;
				return;
			} // if
			if(!is_paused_) {
				stoppoint_ = boost::chrono::steady_clock::now();
				chstop_ = stoppoint_.time_since_epoch();
				stop_ = chstop_.count();
			} // if
			boost::chrono::duration<long long, boost::nano> temp1 = stoppoint_ - startpoint_;
			elapsed_ = temp1.count();
			is_running_ = false;
		} // stop()

		double lap() {
			if(!is_running_) {
				std::cerr << "error: timer is not running" << std::endl;
				return 0.0;
			} // if
			boost::chrono::steady_clock::time_point temppoint = boost::chrono::steady_clock::now();
			boost::chrono::duration<long long, boost::nano> temp1 = temppoint - startpoint_;
			startpoint_ = temppoint;
			elapsed_ += temp1.count();
			return temp1.count();
		} // lap()

		void pause() {
			if(!is_running_) {
				std::cerr << "error: timer is not running" << std::endl;
				return;
			} // if
			if(is_paused_) {
				std::cerr << "warning: timer is already paused. ignoring pause" << std::endl;
				return;
			} // if
			stoppoint_ = boost::chrono::steady_clock::now();
			chstop_ = stoppoint_.time_since_epoch();
			stop_ = chstop_.count();
			is_paused_ = true;
		} // pause()

		void resume() {
			if(!is_running_) {
				std::cerr << "error: timer is not running" << std::endl;
				return;
			} // if
			if(!is_paused_) {
				std::cerr << "warning: timer is not paused. ignoring resume" << std::endl;
				return;
			} // if
			boost::chrono::steady_clock::time_point temppoint = boost::chrono::steady_clock::now();
			boost::chrono::duration<long long, boost::nano> temp1 = temppoint - stoppoint_;
			startpoint_ += temp1;
			chstart_ = startpoint_.time_since_epoch();
			start_ = chstart_.count();
			is_paused_ = false;
		} // resume()


		double elapsed_sec() { return elapsed_ / 1e9; }
		double elapsed_msec() { return elapsed_ / 1e6; }
		double elapsed_usec() { return elapsed_ / 1e3; }
		double elapsed_nsec() { return elapsed_; }

}; // class BoostChronoTimer

} // namespace woo


#endif // _WOOBOOSTCHRONOTIMERS_HPP_
