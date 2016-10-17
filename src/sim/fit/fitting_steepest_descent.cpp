/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: fitting_steepest_descent.cpp
 *  Created: Dec 06, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

//#include <gsl/gsl_linalg.h>

#include "hipgisaxs_main.hpp"


namespace hig {

  // for now m = n = 2
  bool mldivide(int m, real_t* a_mat, real_t* &c_mat) {
  /*  double *a_data = new (std::nothrow) double[m * m];
    for(int i = 0; i < m; ++ i) {
      for(int j = 0; j < m; ++ j) {
        a_data[m * i + j] = (double) a_mat[m * i + j];
      } // for
    } // for
    gsl_matrix_view a = gsl_matrix_view_array(a_data, m, m);
                    
    int s;
    gsl_matrix *x = gsl_matrix_alloc(m, m);
    gsl_permutation* p = gsl_permutation_alloc(m);
    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_invert(&a.matrix, p, x);

    c_mat = new (std::nothrow) real_t[m * m];
    for(int i = 0; i < m; ++ i) {
      for(int j = 0; j < m; ++ j) {
        c_mat[m * i + j] = (real_t) gsl_matrix_get(x, i, j);
      } // for
    } // for
                                     
    gsl_permutation_free(p);
    gsl_matrix_free(x);
    */

    // shortcut for 2x2 matrix a_mat
    if(m != 2) {
      std::cerr << "error: only 2x2 matrix case is implemented" << std::endl;
      return false;
    } // if
    real_t det_a = a_mat[0] * a_mat[3] - a_mat[1] * a_mat[2];
    c_mat = new (std::nothrow) real_t[m * m];
    c_mat[0] = a_mat[3] / det_a;
    c_mat[1] = - a_mat[1] / det_a;
    c_mat[2] = - a_mat[2] / det_a;
    c_mat[3] = a_mat[0] / det_a;

    return true;
  } // mldivide()

  // currently only hacked for spheres, with radius and sd as two parameters
  bool HipGISAXS::fit_steepest_descent(real_t zcut,
          real_t radius_min, real_t radius_max, real_t radius_num,
          real_t sd_min, real_t sd_max, real_t sd_num,
          unsigned int dim, MPI::Intracomm& world_comm,
          int x_min, int x_max, int x_step) {
    int mpi_rank = world_comm.Get_rank();

    if(!init_steepest_fit(world_comm, zcut)) return false;

    int num_alphai = 0, num_phi = 0, num_tilt = 0;;

    real_t alphai_min, alphai_max, alphai_step;
    HiGInput::instance().scattering_alphai(alphai_min, alphai_max, alphai_step);
    if(alphai_max < alphai_min) alphai_max = alphai_min;
    if(alphai_min == alphai_max || alphai_step == 0) num_alphai = 1;
    else num_alphai = (alphai_max - alphai_min) / alphai_step + 1;

    real_t phi_min, phi_max, phi_step;
    HiGInput::instance().scattering_inplanerot(phi_min, phi_max, phi_step);
    if(phi_step == 0) num_phi = 1;
    else num_phi = (phi_max - phi_min) / phi_step + 1;

    real_t tilt_min, tilt_max, tilt_step;
    HiGInput::instance().scattering_tilt(tilt_min, tilt_max, tilt_step);
    if(tilt_step == 0) num_tilt = 1;
    else num_tilt = (tilt_max - tilt_min) / tilt_step + 1;

    std::cout << "**                    Num alphai: " << num_alphai << std::endl
          << "**                       Num phi: " << num_phi << std::endl
          << "**                      Num tilt: " << num_tilt << std::endl;

    // prepare parameters

    std::vector<std::vector<real_t> > params;
    int num_params = 2;
    std::vector<real_t> temp;
    real_t deltap = 0.0;
    if(radius_num <= 1)
      temp.push_back(radius_min);
    else {
      deltap = fabs(radius_max - radius_min) / (radius_num - 1);
      for(int i = 0; i < radius_num; ++ i) {
        temp.push_back(radius_min + i * deltap);
      } // for
    } // if-else
    params.push_back(temp);
    temp.clear();
    if(sd_num <= 1)
      temp.push_back(sd_min);
    else {
      deltap = fabs(sd_max - sd_min) / (sd_num - 1);
      for(int i = 0; i < sd_num; ++ i) {
        temp.push_back(sd_min + i * deltap);
      } // for
    } // if-else
    params.push_back(temp);
    temp.clear();

    // this will work only on one shape and one structure

    const real_t err_threshold = 1e-8;
    const unsigned int max_iter = 200;

    std::vector<real_t> param_vals;
    //param_vals.push_back(16.0);
    //param_vals.push_back(6.0);
    param_vals.push_back(23.0);
    param_vals.push_back(2.0);
    std::vector<real_t> param_deltas;
    param_deltas.push_back(0.05);
    param_deltas.push_back(0.05);
    real_t gamma_const = 0.05;

    real_t qdeltay = QGrid::instance().delta_y();

    real_t alpha_i = alphai_min;
    // high level of parallelism here (alphai, phi, tilt) for dynamicity ...
    for(int i = 0; i < num_alphai; i ++, alpha_i += alphai_step) {
      real_t alphai = alpha_i * PI_ / 180;
      real_t phi = phi_min;
      for(int j = 0; j < num_phi; j ++, phi += phi_step) {
        real_t tilt = tilt_min;
        for(int k = 0; k < num_tilt; k ++, tilt += tilt_step) {

          std::cout << "-- Computing reference GISAXS "
                << i * num_phi * num_tilt + j * num_tilt + k + 1 << " / "
                << num_alphai * num_phi * num_tilt
                << " [alphai = " << alpha_i << ", phi = " << phi
                << ", tilt = " << tilt << "] ..." << std::endl;

          /* run the reference gisaxs simulation using input params */
          real_t* ref_data = NULL;
          if(!run_gisaxs(alpha_i, alphai, phi, tilt, ref_data, world_comm)) {
            if(mpi_rank == 0) std::cerr << "error: could not finish successfully" << std::endl;
            return false;
          } // if

          if(dim != 1) {
            std::cerr << "uh-oh: only 1D is supported for now" << std::endl;
            return false;
          } // if

          real_t* ref_z_cut = new (std::nothrow) real_t[nqy_];
          for(unsigned int iy = 0; iy < nqy_; ++ iy) {
            // assuming nqz_ == 1 ...
            ref_z_cut[iy] = ref_data[nqx_ * iy + 0];
          } // for

          delete[] ref_data;

          // this will store z cut values for each iteration for plotting later
          real_t* z_cuts = new (std::nothrow) real_t[nqy_ * max_iter];
          real_t* temp_zcuts = new (std::nothrow) real_t[nqy_];

          // do some preprocessing
          // start the main loop, bound by max_iter and err_threshold
          //   compute gisaxs for current parameter values
          //   compute the neighbors parameter values
          //   for 12 combinations of current and neighbors, compute gisaxs and error
          //   compute the derivatives (gradient) and error stuff
          //   update parameter values
          // compute the error surface

          real_t err = 10.0;
          std::vector<real_t> param1_list;
          std::vector<real_t> param2_list;
          structure_iterator_t structure_iter = HiGInput::instance().structure_begin();
          Structure* structure = &((*structure_iter).second);
          Shape* shape = HiGInput::instance().shape(*structure);
          shape_param_iterator_t shape_param = (*shape).param_begin();
          real_t* data = NULL;
          std::vector<real_t> param_error_data;
          for(unsigned int iter = 0; iter < max_iter; ++ iter) {
            param1_list.clear();
            param1_list.push_back(param_vals[0] - 2 * param_deltas[0]);  // p1mm
            param1_list.push_back(param_vals[0] - param_deltas[0]);    // p1m
            param1_list.push_back(param_vals[0]);            // p1
            param1_list.push_back(param_vals[0] + param_deltas[0]);    // p1p
            param1_list.push_back(param_vals[0] + 2 * param_deltas[0]);  // p1pp
            param2_list.clear();
            param2_list.push_back(param_vals[1] - 2 * param_deltas[1]);  // p2mm
            param2_list.push_back(param_vals[1] - param_deltas[1]);    // p2m
            param2_list.push_back(param_vals[1]);            // p2
            param2_list.push_back(param_vals[1] + param_deltas[1]);    // p2p
            param2_list.push_back(param_vals[1] + 2 * param_deltas[1]);  // p2pp

            // current point
            (*shape_param).second.mean(param1_list[2]);
            (*shape_param).second.deviation(param2_list[2]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              z_cuts[iter * nqy_ + iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err22 = compute_cut_fit_error(z_cuts + iter * nqy_, ref_z_cut, qdeltay);

            // 12 neighbors

            (*shape_param).second.mean(param1_list[0]);
            (*shape_param).second.deviation(param2_list[2]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err02 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[1]);
            (*shape_param).second.deviation(param2_list[1]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err11 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[1]);
            (*shape_param).second.deviation(param2_list[2]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err12 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[1]);
            (*shape_param).second.deviation(param2_list[3]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err13 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[2]);
            (*shape_param).second.deviation(param2_list[0]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err20 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[2]);
            (*shape_param).second.deviation(param2_list[1]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err21 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[2]);
            (*shape_param).second.deviation(param2_list[3]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err23 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[2]);
            (*shape_param).second.deviation(param2_list[4]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err24 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[3]);
            (*shape_param).second.deviation(param2_list[1]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err31 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[3]);
            (*shape_param).second.deviation(param2_list[2]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err32 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[3]);
            (*shape_param).second.deviation(param2_list[3]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err33 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            (*shape_param).second.mean(param1_list[4]);
            (*shape_param).second.deviation(param2_list[2]);
            if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
              if(mpi_rank == 0)
                std::cerr << "error: could not finish successfully" << std::endl;
              return false;
            } // if
            for(unsigned int iy = 0; iy < nqy_; ++ iy) {
              // assuming nqz_ == 1 ...
              temp_zcuts[iy] = data[nqx_ * iy];
            } // for
            delete[] data; data = NULL;
            real_t err42 = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);

            // 22  0
            // 02  1mm
            // 11  1m2m
            // 12  1m
            // 13  1m2p
            // 20  2mm
            // 21  2m
            // 23  2p
            // 24  2pp
            // 31  1p2m
            // 32  1p
            // 33  1p2p
            // 42  1pp

            real_t derr1 = (err32 - err12) / (2 * param_deltas[0]);
            real_t derr2 = (err23 - err21) / (2 * param_deltas[1]);
            err = sqrt(derr1 * derr1 + derr2 * derr2);
            std::cout << "++ Iteration: " << iter << ", Error: " << err << std::endl;
            std::cout << "++ Parameter 1: " << param_vals[0]
                  << ", Parameter 2: " << param_vals[1] << std::endl;
            param_error_data.push_back(iter);
            param_error_data.push_back(param_vals[0]);
            param_error_data.push_back(param_vals[1]);
            param_error_data.push_back(err);
            if(err < err_threshold) break;

            real_t herr11 = (err42 + err02 - 2 * err22) /
                      (4 * param_deltas[0] * param_deltas[0]);
            real_t herr12 = (err33 - err13 - (err31 - err11)) /
                      (4 * param_deltas[0] * param_deltas[1]);
            real_t herr21 = (err33 - err13 - (err31 - err11)) /
                      (4 * param_deltas[0] * param_deltas[1]);
            real_t herr22 = (err24 + err20 - 2 * err22) /
                      (4 * param_deltas[1] * param_deltas[1]);
            real_t* herr = new (std::nothrow) real_t[2 * 2];
            herr[0] = herr11;
            herr[1] = herr12;
            herr[2] = herr21;
            herr[3] = herr22;
            real_t* herrinv;
            mldivide(2, herr, herrinv);

            param_vals[0] -= gamma_const * (herrinv[0] * derr1 + herrinv[1] * derr2);
            param_vals[1] -= gamma_const * (herrinv[2] * derr1 + herrinv[3] * derr2);

            delete[] herrinv;
            delete[] herr;
          } // for

          // compute the error surface
          std::vector<std::vector<real_t> >::iterator mean_iter = params.begin();
          std::vector<std::vector<real_t> >::iterator sd_iter = mean_iter + 1;
          std::vector<real_t> err_surface;
          for(std::vector<real_t>::iterator curr_mean = (*mean_iter).begin();
              curr_mean != (*mean_iter).end(); ++ curr_mean) {
            for(std::vector<real_t>::iterator curr_sd = (*sd_iter).begin();
                curr_sd != (*sd_iter).end(); ++ curr_sd) {
              (*shape_param).second.mean(*curr_mean);
              (*shape_param).second.deviation(*curr_sd);
              if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
                if(mpi_rank == 0)
                  std::cerr << "error: could not finish successfully" << std::endl;
                return false;
              } // if
              for(unsigned int iy = 0; iy < nqy_; ++ iy) {
                // assuming nqz_ == 1 ...
                temp_zcuts[iy] = data[nqx_ * iy];
              } // for
              delete[] data; data = NULL;
              real_t curr_err = compute_cut_fit_error(temp_zcuts, ref_z_cut, qdeltay);
              err_surface.push_back(*curr_mean);
              err_surface.push_back(*curr_sd);
              err_surface.push_back(curr_err);
            } // for
          } // for

          // write data to files
          // define output filename
          std::stringstream alphai_b, phi_b, tilt_b;
          std::string alphai_s, phi_s, tilt_s;
          alphai_b << alpha_i; alphai_s = alphai_b.str();
          phi_b << phi; phi_s = phi_b.str();
          tilt_b << tilt; tilt_s = tilt_b.str();
          std::string param_error_file(HiGInput::instance().param_pathprefix() +
                        "/" + HiGInput::instance().runname() +
                        "/param_error_ai=" + alphai_s + "_rot=" + phi_s +
                        "_tilt=" + tilt_s + ".dat");
          std::string z_cut_file(HiGInput::instance().param_pathprefix() +
                        "/" + HiGInput::instance().runname() +
                        "/z_cut_ai=" + alphai_s + "_rot=" + phi_s +
                        "_tilt=" + tilt_s + ".dat");
          std::string err_surf_file(HiGInput::instance().param_pathprefix() +
                        "/" + HiGInput::instance().runname() +
                        "/err_surf_ai=" + alphai_s + "_rot=" + phi_s +
                        "_tilt=" + tilt_s + ".dat");
          // write param_error_data
          std::ofstream param_error_f(param_error_file.c_str());
          for(std::vector<real_t>::iterator pei = param_error_data.begin();
              pei != param_error_data.end(); pei += 4) {
            param_error_f << *pei << "\t" << *(pei + 1) << "\t" << *(pei + 2) << "\t"
                    << *(pei + 3) << std::endl;
          } // for
          param_error_f.close();
          // write ref_z_cut and z_cuts
          std::ofstream zcut_f(z_cut_file.c_str());
          for(unsigned int yy = 0; yy < nqy_; ++ yy) {
            zcut_f << ref_z_cut[yy] << "\t";
          } // for
          zcut_f << std::endl;
          for(unsigned int i = 0; i < max_iter; ++ i) {
            for(unsigned int yy = 0; yy < nqy_; ++ yy) {
              zcut_f << z_cuts[i * nqy_ + yy] << "\t";
            } // for
            zcut_f << std::endl;
          } // for
          zcut_f.close();
          // write error surface
          std::ofstream err_surf_f(err_surf_file.c_str());
          for(std::vector<real_t>::iterator surfi = err_surface.begin();
              surfi != err_surface.end(); surfi += 3) {
            err_surf_f << *surfi << "\t" << *(surfi + 1) << "\t" << *(surfi + 2) << std::endl;
          } // for
          err_surf_f.close();

          (*shape_param).second.mean(22.0);
          (*shape_param).second.deviation(7.0);

          param_error_data.clear();
          delete[] temp_zcuts;
          delete[] z_cuts;
          delete[] ref_z_cut;

          std::cout << "parameter values: " << param_vals[0] << ", " << param_vals[1]
                << " [error: " << err << "]" << std::endl;

          // synchronize all procs after each run
          world_comm.Barrier();
        } // for tilt
      } // for phi
    } // for alphai

    return true;
  } // HipGISAXS::fit_all_gisaxs()


  real_t HipGISAXS::compute_cut_fit_error(real_t* cut, real_t* ref, real_t dy) {
    real_t err_sum = 0.0;
    for(unsigned int i = 0; i < nqy_; ++ i) {
      err_sum += fabs(cut[i] - ref[i]);
    } // for
    return err_sum * dy;
  } // HipGISAXS::compute_cut_fit_error()
} // namespace hig
