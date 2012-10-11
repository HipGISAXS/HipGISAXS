/***
  *  $Id: ff.hpp 46 2012-08-23 02:01:21Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff.hpp
  *  Created: Jul 18, 2012
  *  Modified: Tue 09 Oct 2012 12:02:35 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _FF_HPP_
#define _FF_HPP_

#include <string>

#include "typedefs.hpp"
#include "ff_ana.hpp"
#include "ff_num.hpp"


namespace hig {

	// fix the numeric and non-numeric version ... unify them somehow ...
	class FormFactor {
		private:
			bool is_analytic_;
			std::vector <complex_t> ff_;						/* the form factor data */
			AnalyticFormFactor analytic_ff_;

			cucomplex_t* nff_;									/* the cuda version for numeric */
			NumericFormFactor <float_t, cucomplex_t> numeric_ff_; // same as MFormFactor
			//NumericFormFactor <float, cuFloatComplex> numeric_ff_; // same as MFormFactor

		public:
			FormFactor(int a, int b, int c): nff_(NULL), numeric_ff_(a, b, c), is_analytic_(false) { }
			FormFactor(int s): nff_(NULL), numeric_ff_(s), is_analytic_(false) { }
#ifdef KERNEL2
			FormFactor(): nff_(NULL), numeric_ff_(2, 4, 4), is_analytic_(false) { } // default cuda block size
#else
			FormFactor(): nff_(NULL), numeric_ff_(4), is_analytic_(false) { } // default cuda block size
#endif // KERNEL2
			~FormFactor() {
				if(nff_ != NULL) delete[] nff_;
				nff_ = NULL;
			} // ~FormFactor()

			void clear(void);

			bool compute_form_factor(ShapeName shape, std::string shape_filename,
									shape_param_list_t& params,
									float_t single_thickness,
									vector3_t& transvec, float_t shp_tau, float_t shp_eta,
									vector3_t& rot1, vector3_t& rot2, vector3_t& rot3,
									MPI::Intracomm& world_comm);

			//cucomplex_t operator[](unsigned int i) const { return nff_[i]; }
			complex_t operator[](unsigned int i) const {
				if(is_analytic_) return ff_[i];
				else return complex_t(nff_[i].x, nff_[i].y);
			} // operator[]

			// for testing only ... remove ...
			cucomplex_t* nff() { return nff_; }
			bool read_form_factor(const char* filename, unsigned int nqx, unsigned int nqy, unsigned int nqz);
			void print_ff(unsigned int nqx, unsigned int nqy, unsigned int nqz);
			void save_ff(unsigned int nqx, unsigned int nqy, unsigned int nqz, const char* filename);
			void printff(unsigned int nqx, unsigned int nqy, unsigned int nqz) {
				std::cout << "ff:" << std::endl;
				for(unsigned int i = 0; i < nqx * nqy * nqz; ++ i) {
					std::cout << nff_[i].x << "," << nff_[i].y << "\t";
				} // for
				std::cout << std::endl;
			} // printff()
	}; // class FormFactor

} // namespace


#endif // _FF_HPP_
