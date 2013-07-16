/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: layer.cpp
 *  Created: Jun 13, 2012
 *  Modified: Tue 16 Jul 2013 11:51:11 AM PDT
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

#include <iostream>
#include <iomanip>

#include "layer.hpp"

namespace hig {

	Layer::Layer() { }
	Layer::~Layer() { }

	void Layer::init() {
		clear();
	} // Layer::init()

	void Layer::clear() {
		key_.clear();
		order_ = 0;
		thickness_ = 0.0;
		refindex_.delta(0.0); refindex_.beta(0.0);
	} // Layer::clear()

	void Layer::print() {
		std::cout << " key_ = " << key_ << std::endl
					<< " order_ = " << order_ << std::endl
					<< " thickness_ = " << thickness_ << std::endl
					<< " refindex_ = [" << refindex_.delta() << ", "
					<< refindex_.beta() << "]" << std::endl << std::endl;
	} // Layer::print()

} // namespace hig
