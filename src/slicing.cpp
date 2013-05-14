/***
  *  Project:
  *
  *  File: slicing.cpp
  *  Created: May 13, 2013
  *  Modified: Tue 14 May 2013 02:48:38 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "slicing.hpp"

namespace hig {

	bool HiGInput::construct_slices() {

		if(param_nslices() == -1) {	// this case will never happen :-/ ... check ...

			std::cout << "----> SLICING - this is not there in " << __FILE__ << ":" << __LINE__ << std::endl;
			// ...
			return false;
		} else if(param_nslices() < 2) {
			return true;
		} else {
			float_t z_min = min(layer_z_val_min(), domain_.min_[2]);
			float_t z_max = 0.0;		// FIXME ??? need to hard code?

			float_t d_x = fabs(domain_.max_[0] - domain_.min_[0]);
			float_t d_y = fabs(domain_.max_[1] - domain_.min_[1]);
			float_t d_z = fabs(z_max - z_min);

			// slice thickness (uniform slicing)
			float_t d_slice = d_z / (param_nslices() + 1);

			slices_.clear();
			Slice s_vacuum; s_vacuum.clear();
			slices_.push_back(s_vacuum);

			for(int s_i = 0; s_i < param_nslices(); ++ s_i) {
				Slice s;
				s.z_top(z_max + s_i * d_slice);
				s.z_bottom(z_top - d_slice);
				s.volume(d_x * d_y * d_slice);
				if(!s.compute_slice_data(s_i)) return false;
				slices_.push_back(s);
			} // for
		} // if-else

		return true;
	} // HiGInput::construct_slices()


	Slice::Slice() { clear(); }
	Slice::~Slice() { }


	void Slice::init() { clear(); }


	void Slice::clear() {
		z_top_ = 0.0;
		z_bottom_ = 0.0;
		thickness_ = 0.0;
		avg_refindex.delta_ = 0.0;
		avg_refindex.beta_ = 0.0;
	} // Slice::clear()


	bool Slice::compute_slice_data(int s_i, std::vector <float_t> & volumes) {
		// index "s_i" not used at all ??????? .... do we need it???

		// layers stuff
		l_refindex.clear();
		int l_order = HiGInput::instance().num_layers() + 1;
		if(z_bottom_ < 0) {
			l_refindex = HiGInput::instance().substrate_refindex();	// default
			// find the layer contaning z_bottom_
			// TODO: for now assuming that slice is contained in just one layer ...
			for(layer_citerator_t l_i = HiGInput::instance().layers_begin();
					l_i < HiGInput::instance().layers_end(); ++ l_i) {
				if((*l_i).z_val() <= z_bottom_) {
					l_order = (*l_i).order();
					l_refindex = (*l_i).refindex();
					break;
				} // if
			} // for
		} // if

		float_t sample_distr = 0;
		avg_refindex.clear();
		std::vector <float_t> volumes; volumes.clear();

		// structures stuff
		for(structure_iterator_t st = HiGInput::instance().structure_begin();
				st != HiGInput::instance().structure_end(); ++ st) {
			Structure *curr_struct = &((*st).second);
			float_t abs_z = (*curr_struct).absolute_z_val();
			Shape *curr_shape = HiGInput::instance().shape(*curr_struct);

			// shape stuff
			vector3_t struct_min(0, 0, 0), struct_max(0, 0, 0);
			if(!HiGInput::instance().compute_shape_domain(*curr_shape, struct_min, struct_max, 0))
				return false;
				// FIXME: mpi_rank above ...
			vector3_t t = (*curr_struct).grain_transvec();
			struct_min = struct_min + t; struct_max = struct_max + t;
			vector3_t struct_size = (struct_max - struct_min).fabs();
			struct_size[2] = min(struct_size[2], fabs(z_top_ - z_bottom_));
			float_t sample_vol = struct_size.reduce_mul();
			// TODO: handle the case for shapes with triangular base ...
			if((*curr_shape).type() == shape_something || (*curr_shape).type() == shape_anotherthing) {
				// TODO ...
				std::cout << "----> implement me please in " << __FILE__ << ":" << __LINE__ << std::endl;
				return false;
			} // if

			// lattice stuff
			vector3_t grain_repets = (*curr_struct).grain_repetition();
			float_t az = (*curr_struct).lattice_vec_a()[2];
			float_t bz = (*curr_struct).lattice_vec_b()[2];
			float_t cz = (*curr_struct).lattice_vec_c()[2];
			RefractiveIndex curr_refindex = (*curr_struct).grain_refindex();

			float_t struct_vol = 0;

			for(int k = 0; k < grain_repets[2]; ++ k) {
				for(int j = 0; j < grain_repets[1]; ++ j) {
					for(int i = 0; i < grain_repets[0]; ++ i) {
						float_t tz = (i - 1) * az + (j - 1) * bz + (k - 1) * cz + abs_z;
						float_t tz_min = tz + shape_min[2];
						float_t tz_max = tz + shape_max[2];
						if(tz < z_top_ && tz > z_bottom_) {
							avg_refindex = avg_refindex + curr_refindex * sample_vol;
							sample_distr = sample_distr + sample_vol;
							struct_vol = struct_vol + sample_vol;
						} // if
					} // for
				} // for j
			} // for k
			volumes.push_back(struct_vol);
		} // for

		float_t matrix_vol = max(0, slice_vol_ - sample_distr);
		avg_refindex_ = (avg_refindex_ + l_refindex * matrix_vol) / slice_vol_;

		layer_refindex_ = l_refindex;
		sample_volume_ = 0;
		for(std::vector <float_t>::iterator i = volumes.begin(); i != volumes.end(); ++ i)
			sample_volume_ += *i;

		return true;
	} // Slice::compute_slice_data()

} // namespace hig
