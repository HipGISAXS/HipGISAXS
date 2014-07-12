/**
 *	Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *	File: Data.hpp
 *	Created: Dec 26, 2013
 *	Modified: Fri 11 Jul 2014 09:07:59 AM PDT
 *	Description: Abstract class for gisaxs data
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

#ifndef _DATA_HPP_
#define _DATA_HPP_

#include <common/typedefs.hpp>

namespace hig{

	/* image point coords, parallel and vertical */
	typedef	struct {
		float_t crd_par;
		float_t crd_ver;
	}	crd2_t;

	/* image point coord indices, parallel and vertical */
	typedef	struct {
		float_t crd_ind_par;
		float_t crd_ind_ver;
	}	crd2_ind_t;

	/* frame rectangle values & indices to select a subset of image	*/
	typedef struct{
		crd2_t min_pt_;
		crd2_ind_t min_pt_ind_;
		crd2_t max_pt_;
		crd2_ind_t max_pt_ind_;
	} frame_t;

	class Data{

	protected :
		int key_;
		string_t filename_;
		int dim_;
	public:
		Data(){}
		~Data(){}
		virtual bool init(){return false;}

		string_t get_name() const { return filename_;}

		float_vec_t create_vec(float v0, float vN, float& stp, int& nv );

	}; /* class Data */

} /* namespace hig */

#endif /* DATA_HPP_ */
