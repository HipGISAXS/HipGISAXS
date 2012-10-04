/***
  *  $Id: layer.cpp 27 2012-07-15 05:37:29Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: layer.cpp
  *  Created: Jun 13, 2012
  *  Modified: Mon 01 Oct 2012 11:14:49 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
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
