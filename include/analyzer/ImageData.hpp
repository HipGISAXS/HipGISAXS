/**
 *	Project: HipGISAXS
 *
 *	File: ImageData.hpp
 *	Created: December 26, 2013
 *
 *	Author: Slim Chourou <stchourou@lbl.gov>
 *	Developers: Slim Chourou <stchourou@lbl.gov>
 *				Abhinav Sarje <asarje@lbl.gov>
 *				Alexander Hexemer <ahexemer@lbl.gov>
 *				Xiaoye Li <xsli@lbl.gov>
 *
 *	Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *	used by employees of academic research institutions, not-for-profit
 *	research laboratories, or governmental research facilities. Please read the
 *	accompanying LICENSE file before downloading the software. By downloading
 *	the software, you are agreeing to be bound by the terms of this
 *	NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _IMAGEDATA_HPP_
#define _IMAGEDATA_HPP_

#include <iostream>

#include <common/enums.hpp>
#include <analyzer/enums.hpp>
#include <analyzer/Data.hpp>

namespace hig{

	class ImageData : public Data {

	private :
		OutputRegionType mode_;
		float_vec_t data_;
		float_vec_t axis_par_;
		float_vec_t axis_ver_;
		unsigned int n_par_;
		unsigned int n_ver_;

		float_vec_t read_string_values(string_t line);

	public:
		ImageData() { }

		ImageData(std::string filename) {
			if(!read(filename)) exit(2);
		} // ImageData()

		ImageData(const float_vec_t& data, float_vec_t axis_par, float_vec_t axis_ver,
				OutputRegionType type) :
				data_(data), axis_par_(axis_par), axis_ver_(axis_ver), mode_(type) {
			n_par_ = axis_par_.size();
			n_ver_ = axis_ver_.size();
		} // ImageData()

		ImageData(int npar, int nver) {
			n_par_ = npar;
			n_ver_ = nver;
			data_.clear();
		} // ImageData()

		~ImageData() { }

		//virtual bool load(string_t filename){return false;}
		//virtual bool init(){return false;}

		/* setters */
		void set_data(const float_vec_t& data) { data_ = data; }
		void set_data(const float* data) {
			data_.clear();
			for(unsigned int i = 0; i < n_par_ * n_ver_; ++ i) data_.push_back(data[i]);
		} // set_data()

		/* getters */
		//float img_p (int iv, int ip) const;
		//float img_q	(float qv, float qp) const;
		float_t operator()(int i, int j) const;
		float_t operator()(unsigned int i, unsigned int j) const;
		float_t operator()(float_t qi, float_t qj) const;

		//float_mat_t img() const {return img_;}
		unsigned int n_par() const { return n_par_;}
		unsigned int n_ver() const { return n_ver_;}

		bool read(string_t filename);
		void save(string_t filename) const;

		//float_t* convert_data();
		float_t* data() { return &data_[0]; }	// TODO: remove this eventually ...
		//unsigned int size() { return n_par_ * n_ver_; }
		unsigned int size() { return data_.size(); }

		void print() const;
	}; /* class ImageData */

} /* namespace hig */

#endif /* IMAGEDATA_HPP_ */
