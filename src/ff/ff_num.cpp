/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num.cpp
 *  Created: Jul 18, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *          Dinesh Kumar  <dkumar@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <algorithm>
//#if (defined(__SSE3__) || defined(INTEL_SB_AVX)) && !defined(USE_GPU) && !defined(__APPLE__)
//  #include <malloc.h>
//#endif

#include <woo/timer/woo_boostchronotimers.hpp>

#include <ff/ff_num.hpp>
#include <common/parameters.hpp>
#include <common/cpu/parameters_cpu.hpp>
#include <model/qgrid.hpp>
#include <utils/utilities.hpp>
#include <file/objectshape_reader.hpp>
#include <file/rawshape_reader.hpp>

namespace hig {

  bool NumericFormFactor::init(RotMatrix_t& rot, std::vector<complex_t>& ff) {
    nqy_ = QGrid::instance().nqy();
    nqz_ = QGrid::instance().nqz_extended();
    ff.clear();
    rot_ = rot;
    return true;
  } // NumericFormFactor::init()


  // default (new)
  bool NumericFormFactor::compute2(const char * filename, complex_vec_t &ff,
                                   RotMatrix_t & rot
                                   #ifdef USE_MPI
                                    , woo::MultiNode & world_comm, std::string comm_key
                                   #endif
                                   ) {
    // initialize 
    init(rot, ff);

    // read file
    std::vector<vertex_t> vertices;
    std::vector<std::vector<int> > faces;
    std::vector<std::vector<int> > dummy;
    ObjectShapeReader shape_reader;
    if(!shape_reader.load_object(filename, vertices, faces, dummy)) {
        std::cerr << "error: shape reader failed to load triangles" << std::endl;
        return false;
    } // if

    // construct triangles
    int num_triangles = faces.size();
    triangle_t* triangles = new (std::nothrow) triangle_t[num_triangles];
    for(int i = 0; i < num_triangles; i++) {
      triangles[i].v1[0] = vertices[faces[i][0]-1].x;
      triangles[i].v1[1] = vertices[faces[i][0]-1].y;
      triangles[i].v1[2] = vertices[faces[i][0]-1].z;

      triangles[i].v2[0] = vertices[faces[i][1]-1].x;
      triangles[i].v2[1] = vertices[faces[i][1]-1].y;
      triangles[i].v2[2] = vertices[faces[i][1]-1].z;

      triangles[i].v3[0] = vertices[faces[i][2]-1].x;
      triangles[i].v3[1] = vertices[faces[i][2]-1].y;
      triangles[i].v3[2] = vertices[faces[i][2]-1].z;
    } // for

    #ifdef USE_MPI
      int num_procs = world_comm.size(comm_key);
      int rank = world_comm.rank(comm_key);
      bool master = world_comm.is_master(comm_key);
    #else
      bool master = true;
    #endif // USE_MPI

    if(master) {
      std::cout << "-- Numerical form factor computation ..." << std::endl
                << "**        Using input shape file: " << filename << std::endl
                << "**     Number of input triangles: " << num_triangles << std::endl
                << "**  Q-grid resolution (q-points): " << nqy_ << std::endl
                #ifdef USE_MPI
                << "** Number of processes requested: " << num_procs << std::endl
                #endif
                << std::flush;
    } // if

    // copy q-points
    real_t * qx = new (std::nothrow) real_t[nqy_];
    if(qx == NULL) {
        std::cerr << "error: failure in memory allocation for qx" << std::endl;
        return false;
    } // if
    for(int i = 0; i < nqy_; ++ i) qx[i] = QGrid::instance().qx(i);

    real_t * qy = new (std::nothrow) real_t[nqy_];
    if(qy == NULL) {
        std::cerr << "error: failure in memory allocation for qy" << std::endl;
        return false;
    } // if
    for(int i = 0; i < nqy_; ++ i) qy[i] = QGrid::instance().qy(i);

    #ifdef FF_NUM_GPU   // for GPU
    cucomplex_t * qz = new (std::nothrow) cucomplex_t[nqz_];
    if(qz == NULL) {
      std::cerr << "error: failure in memory allocation for qz" << std::endl;
      return false;
    } // if
    for(int i = 0; i < nqz_; ++ i) {
      qz[i].x = QGrid::instance().qz_extended(i).real();
      qz[i].y = QGrid::instance().qz_extended(i).imag();
    } // for
    #else               // for CPU
    complex_t * qz = new (std::nothrow) complex_t[nqz_];
    if(qz == NULL) {
        std::cerr << "error: failure in memory allocation for qz" << std::endl;
        return false;
    } // if
    for(int i = 0; i < nqz_; ++ i) qz[i] = QGrid::instance().qz_extended(i);
    #endif              // FF_NUM_GPU

    real_t compute_time = 0.0;
    #ifdef FF_NUM_GPU
      cucomplex_t * p_ff = NULL;
      if(num_triangles != gff_.compute_exact_triangle(triangles, num_triangles,
                                                      p_ff, nqy_, qx, qy,
                                                      nqz_, qz, rot_, compute_time)) {
        std::cerr << "error: calculation of numerical form-factor failed" << std::endl;
        return false;
      } // if
      for(int i = 0; i < nqz_; ++ i) ff.push_back(complex_t(p_ff[i].x, p_ff[i].y));
      std::cout << "**        FF GPU compute time: " << compute_time << " ms." << std::endl;
    #else             // on CPU
      complex_t * p_ff = new (std::nothrow) complex_t[nqz_];
      if(p_ff == NULL){
        std::cerr << "error: failed to allocate memory of size "
                  << nqz_ * sizeof(complex_t) << std::endl;
        return false;
      } // if
      if(num_triangles != cff_.compute_exact_triangle(triangles, num_triangles, 
                                                      p_ff, nqy_, qx, qy, 
                                                      nqz_, qz, rot_, compute_time)) {
        std::cerr << "error: calculation of numerical form-factor failed" << std::endl;
        return false;
      } // if
      for(int i = 0; i < nqz_; ++ i) ff.push_back(p_ff[i]);
      std::cout << "**        FF CPU compute time: " << compute_time << " ms." << std::endl;
    #endif            // FF_NUM_GPU

    delete[] qx;
    delete[] qy;
    delete[] qz;
    delete[] triangles;
    if(p_ff != NULL) delete[] p_ff;

    return true;
  } // NumericFormFactor::compute2()


  // old
  bool NumericFormFactor::compute(const char * filename, complex_vec_t & ff,
                                  RotMatrix_t & rot
                                  #ifdef USE_MPI
                                    , woo::MultiNode &world_comm, std::string comm_key
                                  #endif
                                  ) {
    real_t comp_time = 0.0;

    // initialize 
    init(rot, ff);

    unsigned int nqy = QGrid::instance().nqy();
    unsigned int nqz = QGrid::instance().nqz_extended();

    // warning: all procs read the shape file!!!!
    // TODO: improve to parallel IO, or one proc reading and sending to all ...
//    #ifndef __SSE3__
      real_vec_t shape_def;
//    #else
//      #ifdef USE_GPU
//        real_vec_t shape_def;
//      #else
//        real_t* shape_def = NULL;
//      #endif
//    #endif
    // use the new file reader instead ...
    unsigned int num_triangles = read_shapes_file(filename, shape_def);
            // TODO ... <--- sadly all procs read this! IMPROVE!!!
  
    #ifdef USE_MPI
    int num_procs = world_comm.size(comm_key);
    int rank = world_comm.rank(comm_key);
    bool master = world_comm.is_master(comm_key);
    #else
    bool master = true;
    #endif

    if(master) {
      std::cout << "-- Numerical form factor computation ..." << std::endl
            << "**        Using input shape file: " << filename << std::endl
            << "**     Number of input triangles: " << num_triangles << std::endl
            << "**  Q-grid resolution (q-points): " <<  nqz << std::endl
            #ifdef USE_MPI
              << "** Number of processes requested: " << num_procs << std::endl
            #endif
            << std::flush;
    } // if
    if(num_triangles < 1) {
      std::cerr << "error: no triangles found in specified definition file" << std::endl;
      return false;
    } // if
  
    // FIXME: this is a yucky temporary fix ... fix properly ...
    real_t* qx = new (std::nothrow) real_t[nqy]();
    real_t* qy = new (std::nothrow) real_t[nqy]();
    #ifdef FF_NUM_GPU
    cucomplex_t* qz = new (std::nothrow) cucomplex_t[nqz]();
    #else
    complex_t* qz = new (std::nothrow) complex_t[nqz]();
    #endif
   
    // create qy_and qz using qgrid instance
    for(unsigned int i = 0; i < nqy; ++ i) qx[i] = QGrid::instance().qx(i);
    for(unsigned int i = 0; i < nqy; ++ i) qy[i] = QGrid::instance().qy(i);
    for(unsigned int i = 0; i < nqz; ++ i) {
    #ifdef FF_NUM_GPU
      qz[i].x = QGrid::instance().qz_extended(i).real();
      qz[i].y = QGrid::instance().qz_extended(i).imag();
    #else
      qz[i] = QGrid::instance().qz_extended(i);
    #endif
    } // for
      
    #ifdef FF_NUM_GPU
    cucomplex_t *p_ff = NULL;
    #else
    complex_t *p_ff = NULL;
    #endif
  
    real_t kernel_time = 0.;
    unsigned int ret_numtriangles = 0;
    #ifdef FF_NUM_GPU  // use GPU
    ret_numtriangles = gff_.compute_approx_triangle(shape_def, 
            p_ff, nqy, qx, qy, nqz, qz, rot_, kernel_time);
    for (int i = 0; i < nqz; i++) ff.push_back(complex_t(p_ff[i].x, p_ff[i].y));
    std::cout << "**        FF GPU compute time: " << kernel_time << " ms." << std::endl;
    #else  // use only CPU
    ret_numtriangles = cff_.compute_approx_triangle(shape_def, 
            p_ff, nqy, qx, qy, nqz, qz, rot_, kernel_time);
    for (int i = 0; i < nqz; i++) ff.push_back(p_ff[i]);
    std::cout << "**        FF CPU compute time: " << kernel_time << " ms." << std::endl;
    #endif

      if(p_ff != NULL) delete[] p_ff;
      delete[] qz;
      delete[] qy;
      delete[] qx;
  }

#ifdef OLD_Q_GRID
  /**
   * main host function
   */
  bool NumericFormFactor::compute(const char* filename, complex_vec_t& ff,
                  vector3_t& rot1, vector3_t& rot2, vector3_t& rot3
                  #ifdef USE_MPI
                    , woo::MultiNode& world_comm, std::string comm_key
                  #endif
                  ) {
    real_t comp_start = 0.0, comp_end = 0.0, comm_start = 0.0, comm_end = 0.0;
    real_t mem_start = 0.0, mem_end = 0.0;
    real_t comp_time = 0.0, comm_time = 0.0, mem_time = 0.0, kernel_time = 0.0, red_time = 0.0;
    real_t total_start = 0.0, total_end = 0.0, total_time = 0.0;

    woo::BoostChronoTimer maintimer, computetimer;
    woo::BoostChronoTimer commtimer, memtimer;

    unsigned int nqx = QGrid::instance().nqx();
    unsigned int nqy = QGrid::instance().nqy();
    unsigned int nqz = QGrid::instance().nqz_extended();

    #ifdef USE_MPI
      bool master = world_comm.is_master(comm_key);
      commtimer.start();
      world_comm.barrier(comm_key);
      commtimer.stop();
      comm_time += commtimer.elapsed_msec();
    #else
      bool master = true;
    #endif

  
    // warning: all procs read the shape file!!!!
    // TODO: improve to parallel IO, or one proc reading and sending to all ...
//    #ifndef __SSE3__
      real_vec_t shape_def;
//    #else
//      #ifdef USE_GPU
//        real_vec_t shape_def;
//      #else
//        real_t* shape_def = NULL;
//      #endif
//    #endif
    // use the new file reader instead ...
    unsigned int num_triangles = read_shapes_file(filename, shape_def);
            // TODO ... <--- sadly all procs read this! IMPROVE!!!
  
    // TODO: temporary ... remove ...
    std::vector<short int> axes(4);      // axes[i] = j
                        // i: x=0 y=1 z=2
                        // j: 0=a 1=b 2=c
    #ifndef AXIS_ROT
      axes[0] = 0; axes[1] = 1; axes[2] = 2;  // default values
    #else
      find_axes_orientation(shape_def, axes);
    #endif

    #ifdef USE_MPI
      int num_procs = world_comm.size(comm_key);
      int rank = world_comm.rank(comm_key);
    #endif

    if(master) {
      std::cout << "-- Numerical form factor computation ..." << std::endl
            << "**        Using input shape file: " << filename << std::endl
            << "**     Number of input triangles: " << num_triangles << std::endl
            << "**  Q-grid resolution (q-points): " << nqx * nqy * nqz << std::endl
                  << "**               NQX x NQY x NQZ: "
            << nqx << " x " << nqy << " x " << nqz << std::endl
            #ifdef USE_MPI
              << "** Number of processes requested: " << num_procs << std::endl
            #endif
            << std::flush;
    } // if
    if(num_triangles < 1) {
      std::cerr << "error: no triangles found in specified definition file" << std::endl;
      return false;
    } // if
  
    #ifdef USE_MPI
      // decompose along y and z directions into blocks
      int p_y = std::floor(sqrt((real_t) num_procs));  // some procs may be idle ...
      int p_z = num_procs / p_y;
    
      int p_nqx = nqx;
      int p_nqy = nqy / p_y + (((rank / p_z) < (int)nqy % p_y) ? 1 : 0);
      int p_nqz = nqz / p_z + (((rank % p_z) < (int)nqz % p_z) ? 1 : 0);

      commtimer.start();

      int idle = 0;
      if(world_comm.rank(comm_key) >= p_y * p_z) idle = 1;
      std::string real_world("ff_num_real_world");
      world_comm.split(real_world, comm_key, idle);

      commtimer.stop();
      comm_time += commtimer.elapsed_msec();
    #else
      int p_y = 1, p_z = 1;
      int p_nqx = nqx;
      int p_nqy = nqy;
      int p_nqz = nqz;
    #endif // USE_MPI
  
    #ifdef FINDBLOCK
      int block_x = 0, block_y = 0, block_z = 0, block_t = 0;
      int block_x_max = 0, block_y_max = 0, block_z_max = 0, block_t_max = 0;
      block_x_max = (nqx < 400) ? nqx : 400;
      block_y_max = (nqy < 400) ? nqy : 400;
      block_z_max = (nqz < 400) ? nqz : 400;
      block_t_max = (num_triangles < 2500) ? num_triangles : 2500;
      block_t = block_t_max;
      for(block_t = block_t_max; block_t > std::min(99, block_t_max - 1); block_t -= 100) {
      for(block_x = block_x_max; block_x > std::min(3, block_x_max - 1); block_x -= 2) {
      for(block_y = block_y_max; block_y > std::min(3, block_y_max - 1); block_y -= 2) {
      for(block_z = block_z_max; block_z > std::min(3, block_z_max - 1); block_z -= 2) {
    #endif
    
    maintimer.start();

    #ifdef USE_MPI
    if(world_comm.rank(comm_key) < p_y * p_z) {    // only the non-idle processors
      bool master = world_comm.is_master(real_world);
      if(master) {
        std::cout << "++  Number of MPI processes used: "
              << world_comm.size(real_world) << std::endl
              << "++                 MPI grid size: 1 x " << p_y << " x " << p_z
              << std::endl << std::flush;
      } // if

      commtimer.start();

      int rank = world_comm.rank(real_world);
      int size = world_comm.size(real_world);

      // create row-wise and column-wise communicators
      int row = rank / p_z, col = rank % p_z;
      world_comm.split("ff_num_row_comm", real_world, row);
      world_comm.split("ff_num_col_comm", real_world, col);

      // perform MPI scan operation to compute y_offset and z_offset

      unsigned int y_offset = 0, z_offset = 0;
      world_comm.scan_sum("ff_num_col_comm", p_nqy, y_offset);
      world_comm.scan_sum("ff_num_row_comm", p_nqz, z_offset);

      commtimer.stop();
      comm_time += commtimer.elapsed_msec();
  
      y_offset -= p_nqy;
      z_offset -= p_nqz;
    #else
      master = true;
      unsigned int y_offset = 0, z_offset = 0;
      int rank = 0;
      int size = 1;
    #endif // USE_MPI

      memtimer.start();

      // FIXME: this is a yucky temporary fix ... fix properly ...
      real_t* qx = new (std::nothrow) real_t[nqx]();
      real_t* qy = new (std::nothrow) real_t[nqy]();
      #ifdef FF_NUM_GPU
        cucomplex_t* qz = new (std::nothrow) cucomplex_t[nqz]();
      #else
        complex_t* qz = new (std::nothrow) complex_t[nqz]();
      #endif
      // create qy_and qz using qgrid instance
      for(unsigned int i = 0; i < nqx; ++ i) {
        qx[i] = QGrid::instance().qx(i);
      } // for
      for(unsigned int i = 0; i < nqy; ++ i) {
        qy[i] = QGrid::instance().qy(i);
      } // for
      for(unsigned int i = 0; i < nqz; ++ i) {
        #ifdef FF_NUM_GPU
          qz[i].x = QGrid::instance().qz_extended(i).real();
          qz[i].y = QGrid::instance().qz_extended(i).imag();
        #else
          qz[i] = QGrid::instance().qz_extended(i);
        #endif
      } // for
      
      #ifdef USE_MPI
        // create p_ff buffers  <----- TODO: IMPROVE for all procs!!!
        real_t *p_qy = NULL;
        p_qy = new (std::nothrow) real_t[p_nqy]();
        if(p_qy == NULL) { return 0; }
        memcpy(p_qy, (void*) (qy + y_offset), p_nqy * sizeof(real_t));
        #ifdef FF_NUM_GPU
          cucomplex_t *p_qz = NULL;
          p_qz = new (std::nothrow) cucomplex_t[p_nqz]();
          if(p_qz == NULL) { delete[] p_qy; return 0; }
          memcpy(p_qz, (void*) (qz + z_offset), p_nqz * sizeof(cucomplex_t));
        #else // TODO: avoid the following ...
          complex_t *p_qz = NULL;
          p_qz = new (std::nothrow) complex_t[p_nqz]();
          if(p_qz == NULL) { delete[] p_qy; return 0; }
          memcpy(p_qz, (void*) (qz + z_offset), p_nqz * sizeof(complex_t));
        #endif  // FF_NUM_GPU
      #else  // no MPI
        real_t *p_qy = qy;
        #ifdef FF_NUM_GPU
          cucomplex_t *p_qz = qz;
        #else
          complex_t *p_qz = qz;
        #endif // FF_NUM_GPU
      #endif  // USE_MPI
  
      memtimer.stop();
      mem_time += memtimer.elapsed_msec();
    
      // compute local

      #ifdef FF_NUM_GPU
        cucomplex_t *p_ff = NULL;
      #else
        complex_t *p_ff = NULL;
      #endif
  
      computetimer.reset();
      computetimer.start();

      unsigned int ret_numtriangles = 0;

      real_t temp_mem_time = 0.0;

      #ifdef FF_NUM_GPU  // use GPU
        #ifdef FF_NUM_GPU_FUSED
          ret_numtriangles = gff_.compute_form_factor_kb_fused(rank, shape_def, axes, p_ff,
                        qx, p_nqx, p_qy, p_nqy, p_qz, p_nqz, 3,
                        rot_,
                        kernel_time, red_time, temp_mem_time
                        #ifdef FINDBLOCK
                          , block_x, block_y, block_z, block_t
                        #endif
                        );
        #else
          ret_numtriangles = gff_.compute_form_factor_db(rank, shape_def, axes, p_ff,
                        qx, p_nqx, p_qy, p_nqy, p_qz, p_nqz,
                        rot_,
                        kernel_time, red_time, temp_mem_time
                        #ifdef FINDBLOCK
                          , block_x, block_y, block_z, block_t
                        #endif
                        );
        #endif
/*      #elif defined USE_MIC  // use MIC
        #ifndef FF_NUM_MIC_KB
          ret_numtriangles = mff_.compute_form_factor_db(rank, shape_def, p_ff,
                        qx, p_nqx, p_qy, p_nqy, p_qz, p_nqz,
                        rot_,
                        kernel_time, red_time, temp_mem_time
                        #ifdef FINDBLOCK
                          , block_x, block_y, block_z, block_t
                        #endif
                        );
        #else
          ret_numtriangles = mff_.compute_form_factor_kb(rank, shape_def,
                        num_triangles,
                        p_ff,
                        qx, p_nqx, p_qy, p_nqy, p_qz, p_nqz, 3,
                        rot_,
                        kernel_time, red_time, temp_mem_time
                        #ifdef FINDBLOCK
                          , block_x, block_y, block_z, block_t
                        #endif
                        );
        #endif */
      #else  // use only CPU
        ret_numtriangles = cff_.compute_form_factor(rank, shape_def,
//                        #ifdef __SSE3__
//                          num_triangles,
//                        #endif
                        p_ff,
                        qx, p_nqx, p_qy, p_nqy, p_qz, p_nqz,
                        rot_,
                        kernel_time, red_time, temp_mem_time
                        #ifdef FINDBLOCK
                          , block_x, block_y, block_z, block_t
                        #endif
                        );
      #endif

      computetimer.stop();
      comp_time += computetimer.elapsed_msec();
      mem_time += (temp_mem_time / 1000);
  
      // gather everything on proc 0
      if(ret_numtriangles > 0) {
        real_t temp_mem_time = 0.0, temp_comm_time = 0.0;
        construct_ff(p_nqx, p_nqy, p_nqz, nqx, nqy, nqz, p_y, p_z, p_ff, ff,
                #ifdef USE_MPI
                  world_comm, real_world,
                #endif
                temp_mem_time, temp_comm_time);
        mem_time += temp_mem_time;
        comm_time += temp_comm_time;
      } // if
  
      /*if(rank == 0) {
        write_slice_to_file(ff, nqx, nqy, nqz, filename, 0, 0);  // x = 0, y = 1, z = 2
                            // only slice along x implemented for now
      } // if*/
  
      #ifdef USE_MPI
        world_comm.barrier(real_world);
      #endif
  
      memtimer.start();
      
      #ifdef FINDBLOCK
        ff.clear();
      #endif
      if(p_ff != NULL) delete[] p_ff;
      #ifdef USE_MPI
        delete[] p_qz;
        delete[] p_qy;
      #endif
      delete[] qz;
      delete[] qy;
      delete[] qx;

      memtimer.stop();
      maintimer.stop();
  
      total_time = maintimer.elapsed_msec();
      mem_time += memtimer.elapsed_msec();
  
      if(master) {
        #ifdef TIME_DETAIL_1
          std::cout
            << "**                FF kernel time: " << kernel_time << " ms." << std::endl
            << "**               FF compute time: " << computetimer.elapsed_msec() << " ms."
            << std::endl
            << "**         FF memory and IO time: " << mem_time * 1000 << " ms." << std::endl
            << "**            Communication time: " << comm_time * 1000 << " ms." << std::endl
            << "**                 Total FF time: " << maintimer.elapsed_msec() << " ms."
            << std::endl << std::flush;
        #endif // TIME_DETAIL_1

        double mflop = 0.0; real_t gflops = 0.0;

        #ifdef USE_GPU
          // flop count for GPU
          //mflop = (double) nqx * nqy * nqz * (42 * num_triangles + 2) / 1000000;
          mflop = (double) nqx * nqy * nqz * (69 * num_triangles + 52) / 1000000;
        #elif defined USE_MIC
          // flop count for MIC
          //mflop = (double) nqx * nqy * nqz * (78 * num_triangles + 18) / 1000000;
          mflop = (double) nqx * nqy * nqz * (111 * num_triangles + 50) / 1000000;
        #elif defined INTEL_SB_AVX
          // flop count for Sandy Bridge with AVX
          // TODO: recount flops ...
          mflop = (double) nqx * nqy * nqz * (85 * num_triangles + 16) / 1000000;
        #else
          // flop count for SSE3 CPU (hopper)
          // TODO: recount flops ...
          mflop = (double) nqx * nqy * nqz * (68 * num_triangles + 20) / 1000000;
        #endif
        //gflops = nidle_num_procs * mflop / kernel_time;
        gflops = mflop / kernel_time;
        std::cout << "**            Kernel performance: " << gflops << " GFLOPS/s" << std::endl;
      } // if
    #ifdef USE_MPI
      world_comm.free("ff_num_row_comm");
      world_comm.free("ff_num_col_comm");
    } // if
    #endif

    #ifdef USE_MPI
      world_comm.barrier(comm_key);
      world_comm.free(real_world);
    #endif

    #ifdef FINDBLOCK
      } // block_t
      } // block_z
      } // block_y
      } // block_x
    #endif

    return true;
  } // NumericFormFactor::compute()

#endif // OLD_Q_GRID

  /**
   * Function to gather partial FF arrays from all processes to construct the final FF.
   * This is a bottleneck for large num procs ...
   */
  bool NumericFormFactor::construct_ff(int p_nqx, int p_nqy, int p_nqz,
                      int nqx, int nqy, int nqz,
                      int p_y, int p_z,
                      #ifdef FF_NUM_GPU
                        cucomplex_t* p_ff,
                      #else
                        complex_t* p_ff,
                      #endif
                      complex_vec_t& ff,
                      #ifdef USE_MPI
                        woo::MultiNode& world_comm, std::string comm_key,
                      #endif
                      real_t& mem_time, real_t& comm_time) {
    real_t mem_start = 0, mem_end = 0, comm_start = 0, comm_end = 0;
    woo::BoostChronoTimer memtimer, commtimer;
    mem_time = 0; comm_time = 0;

    #ifdef USE_MPI
      bool master = world_comm.is_master(comm_key);
      int size = world_comm.size(comm_key);
      int rank = world_comm.rank(comm_key);
    #else
      bool master = true;
      int size = 1;
      int rank = 0;
    #endif

    memtimer.start();
  
    int local_qpoints = p_nqx * p_nqy * p_nqz;
    unsigned long int total_qpoints = nqx * nqy * nqz;
  
    // process 0 creates the main ff, and collects computed p_ff from all others (just use gather)
    ff.clear();
    #ifdef FF_NUM_GPU
      cucomplex_t* all_ff = NULL;    // TODO: improve this ...
    #else
      complex_t* all_ff = NULL;
    #endif
    if(master) {
      ff.reserve(total_qpoints);
      ff.assign(total_qpoints, complex_t(0.0, 0.0));
      #ifdef FF_NUM_GPU
        all_ff = new (std::nothrow) cucomplex_t[total_qpoints];
      #else
        all_ff = new (std::nothrow) complex_t[total_qpoints];
      #endif
    } // if
  
    mem_time += memtimer.elapsed_msec();

    int *recv_p_nqy = new (std::nothrow) int[p_y]();
    recv_p_nqy[0] = p_nqy;
    int *off_p_nqy = new (std::nothrow) int[p_y]();
    off_p_nqy[0] = 0;

    #ifdef FF_NUM_GPU
      cucomplex_t *ff_buffer = new (std::nothrow) cucomplex_t[total_qpoints];
    #else
      complex_t *ff_buffer = new (std::nothrow) complex_t[total_qpoints];
    #endif
    if(ff_buffer == NULL) {
      std::cerr << "error: failed to allocate memory for ff buffer" << std::endl;
      return false;
    } // if

    #ifdef USE_MPI
      // construct stuff for gatherv
      int *recv_counts = new (std::nothrow) int[size]();
      int *displs = new (std::nothrow) int[size]();

      commtimer.start();

      //comm.Allgather(&local_qpoints, 1, MPI::INT, recv_counts, 1, MPI::INT);
      world_comm.allgather(comm_key, &local_qpoints, 1, recv_counts, 1);

      commtimer.stop();
      comm_time += commtimer.elapsed_msec();
      memtimer.start();

      displs[0] = 0;
      for(int i = 1; i < size; ++ i) {
        displs[i] = displs[i - 1] + recv_counts[i - 1];
      } // for
      complex_t *cast_p_ff, *cast_ff;
      #ifdef FF_NUM_GPU
        cast_p_ff = reinterpret_cast<complex_t*>(p_ff);
        cast_ff = reinterpret_cast<complex_t*>(ff_buffer);
      #else
        cast_p_ff = p_ff;
        cast_ff = ff_buffer;
      #endif

      memtimer.stop();
      mem_time += memtimer.elapsed_msec();
  
      commtimer.start();

      world_comm.gatherv(comm_key, cast_p_ff, local_qpoints, cast_ff, recv_counts, displs);
  
      world_comm.gather("ff_num_col_comm", &p_nqy, 1, recv_p_nqy, 1);
    
      commtimer.stop();
      comm_time += commtimer.elapsed_msec();

      for(int i = 1; i < p_y; ++ i) off_p_nqy[i] = off_p_nqy[i - 1] + recv_p_nqy[i - 1];
    #else
      #ifdef FF_NUM_GPU
        memcpy(ff_buffer, p_ff, total_qpoints * sizeof(cucomplex_t));
      #else
        memcpy(ff_buffer, p_ff, total_qpoints * sizeof(complex_t));
      #endif
    #endif // USE_MPI
  
    memtimer.start();

    // move all the data to correct places
    if(rank == 0) {
      unsigned long int ff_index = 0;
      for(int i_nqz = 0; i_nqz < nqz; ++ i_nqz) {
        for(int i_py = 0; i_py < p_y; ++ i_py) {
          unsigned long int ffb_index = nqx * (i_nqz * recv_p_nqy[i_py] +
                              nqz * off_p_nqy[i_py]);
          #ifdef FF_NUM_GPU
            memcpy(&all_ff[ff_index], &ff_buffer[ffb_index],
                nqx * recv_p_nqy[i_py] * sizeof(cucomplex_t));
          #else
            memcpy(&all_ff[ff_index], &ff_buffer[ffb_index],
                nqx * recv_p_nqy[i_py] * sizeof(complex_t));
          #endif
          ff_index += nqx * recv_p_nqy[i_py];
        } // for i_py
      } // for i_nqz
      // put into the final ff buffer
      #ifdef FF_NUM_GPU
        ff.assign(reinterpret_cast<complex_t*>(all_ff),
              reinterpret_cast<complex_t*>(all_ff + total_qpoints));
      #else
        ff.assign(all_ff, all_ff + total_qpoints);
      #endif
    } // if
  
    delete[] ff_buffer;
    #ifdef USE_MPI
      delete[] displs;
      delete[] recv_counts;
    #endif
    delete[] off_p_nqy;
    delete[] recv_p_nqy;
    delete[] all_ff;

    memtimer.stop();
    mem_time += memtimer.elapsed_msec();

    return true;
  } // NumericFormFactor::construct_ff()
  
  
  /**
   * Function to read the input shape file.
   */
  unsigned int NumericFormFactor::read_shapes_file_dat(const char* filename, real_vec_t &shape_def) {
    std::ifstream f(filename);
    if(!f.is_open()) {
      std::cout << "error: could not open file " << filename << std::endl;
      return 0;
    } // if

    real_t s = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;
    while(true) {
      f >> s;
      if(f.eof() || !f.good()) break;
      f >> nx; f >> ny; f >> nz;
      f >> cx; f >> cy; f >> cz;
      shape_def.push_back(s);
      shape_def.push_back(nx);
      shape_def.push_back(ny);
      shape_def.push_back(nz);
      shape_def.push_back(cx);
      shape_def.push_back(cy);
      shape_def.push_back(cz);
    } // while
  
    f.close();
    return shape_def.size() / 7;
  } // NumericFormFactor::read_shapes_file_dat()
  
  
  void NumericFormFactor::find_axes_orientation(real_vec_t &shape_def, std::vector<short int> &axes) {
    real_t min_a = shape_def[4], max_a = shape_def[4];
    real_t min_b = shape_def[5], max_b = shape_def[5];
    real_t min_c = shape_def[6], max_c = shape_def[6];
  
    for(unsigned int i = 0; i + 6 < shape_def.size(); i += 7) {
      min_a = (min_a > shape_def[i + 4]) ? shape_def[i + 4] : min_a ;
      max_a = (max_a < shape_def[i + 4]) ? shape_def[i + 4] : max_a ;
      min_b = (min_b > shape_def[i + 5]) ? shape_def[i + 5] : min_b ;
      max_b = (max_b < shape_def[i + 5]) ? shape_def[i + 5] : max_b ;
      min_c = (min_c > shape_def[i + 6]) ? shape_def[i + 6] : min_c ;
      max_c = (max_c < shape_def[i + 6]) ? shape_def[i + 6] : max_c ;
    } // for
  
    real_t diff_a = max_a - min_a;
    real_t diff_b = max_b - min_b;
    real_t diff_c = max_c - min_c;
  
    // axes[i] = j
    // i: x=0 y=1 z=2
    // j: 0=a 1=b 2=c

    //std::cout << "++ diff_a = " << diff_a << ", diff_b = " << diff_b
    //      << ", diff_c = " << diff_c << std::endl;

    real_vec_t min_point, max_point;
  
    // the smallest one is x, other two are y and z
    if(diff_a < diff_b) {
      if(diff_a < diff_c) {
        // x is a
        axes[0] = 0; axes[1] = 1; axes[2] = 2;
        min_point.push_back(min_a); min_point.push_back(min_b); min_point.push_back(min_c);
        max_point.push_back(max_a); max_point.push_back(max_b); max_point.push_back(max_c);
      } else {
        // x is c
        axes[0] = 2; axes[1] = 0; axes[2] = 1;
        min_point.push_back(min_c); min_point.push_back(min_a); min_point.push_back(min_b);
        max_point.push_back(max_c); max_point.push_back(max_a); max_point.push_back(max_b);
      } // if-else
    } else {
      if(diff_b < diff_c) {
        // x is b
        axes[0] = 1; axes[1] = 0; axes[2] = 2;
        min_point.push_back(min_b); min_point.push_back(min_a); min_point.push_back(min_c);
        max_point.push_back(max_b); max_point.push_back(max_a); max_point.push_back(max_c);
      } else {
        // x is c
        axes[0] = 2; axes[1] = 0; axes[2] = 1;
        min_point.push_back(min_c); min_point.push_back(min_a); min_point.push_back(min_b);
        max_point.push_back(max_c); max_point.push_back(max_a); max_point.push_back(max_b);
      } // if-else
    } // if-else

    std::cout << "++ Shape min point: " << min_point[0] << ", "
          << min_point[1] << ", " << min_point[2] << std::endl;
    std::cout << "++ Shape max point: " << max_point[0] << ", "
          << max_point[1] << ", " << max_point[2] << std::endl;
    std::cout << "++ Shape dimensions: "
          << fabs(max_point[0] - min_point[0]) << " x "
          << fabs(max_point[1] - min_point[1]) << " x "
          << fabs(max_point[2] - min_point[2]) << std::endl;
  } // NumericFormFactor::find_axes_orientation()


  ShapeFileType NumericFormFactor::get_shapes_file_format(const char* filename) {
    std::istringstream file(filename);
    std::string s;    
    while(std::getline(file, s, '.'))// std::cout << s << std::endl;
    if(s.compare("") == 0) return shape_file_null;
    if(s.compare("dat") == 0) return shape_file_data;
    if(s.compare("Dat") == 0) return shape_file_data;
    if(s.compare("DAT") == 0) return shape_file_data;
    if(s.compare("hd5") == 0) return shape_file_hdf5;
    if(s.compare("Hd5") == 0) return shape_file_hdf5;
    if(s.compare("HD5") == 0) return shape_file_hdf5;
    if(s.compare("hdf5") == 0) return shape_file_hdf5;
    if(s.compare("Hdf5") == 0) return shape_file_hdf5;
    if(s.compare("HDF5") == 0) return shape_file_hdf5;
    if(s.compare("obj") == 0) return shape_file_object;
    if(s.compare("Obj") == 0) return shape_file_object;
    if(s.compare("OBJ") == 0) return shape_file_object;
    return shape_file_error;
  } // NumericFormFactor::get_shape_file_format()
  
  
  /**
   * Function to read the shape definition input file in HDF5 format.
   */
  unsigned int NumericFormFactor::read_shapes_file(const char* filename,
//                          #ifndef __SSE3__
                            real_vec_t &shape_def
//                          #else
//                            #ifdef USE_GPU
//                              real_vec_t &shape_def
//                            #else
//                              real_t* &shape_def
//                            #endif
//                          #endif
                          ) {
    unsigned int num_triangles = 0;
    double* temp_shape_def = NULL;
  
    // TODO: shape definition is already in HigInput ...
    // utilize ...
    ShapeFileType type = get_shapes_file_format(filename);
    if(type == shape_file_data) {
      RawShapeReader temp(filename, temp_shape_def, num_triangles);
    } else if(type == shape_file_object) {
      ObjectShapeReader temp(filename, temp_shape_def, num_triangles);
    } else if(type == shape_file_hdf5) {
      #ifdef USE_PARALLEL_HDF5
        h5_shape_reader(filename, &temp_shape_def, &num_triangles);
      #else
        std::cerr << "error: use of parallel hdf5 format has not been enabled in your installation. "
                  << "Please reinstall with the support enabled." << std::endl;
        return false;
      #endif
    } else if(type == shape_file_null) {
      std::cerr << "error: shape definition file extension is null" << std::endl;
      return 0;
    } else if(type == shape_file_error) {
      std::cerr << "error: shape definition file format unknown" << std::endl;
      return 0;
    } else {
      std::cerr << "error: shape definition file format unknown" << std::endl;
      return 0;
    } // if-else

    #ifdef FF_NUM_GPU
      #ifndef KERNEL2
        for(unsigned int i = 0; i < num_triangles * 7; ++ i)
          shape_def.push_back((real_t)temp_shape_def[i]);
      #else // KERNEL2
        for(unsigned int i = 0, j = 0; i < num_triangles * T_PROP_SIZE_; ++ i) {
          if((i + 1) % T_PROP_SIZE_ == 0) shape_def.push_back((real_t) 0.0);  // padding
          else { shape_def.push_back((real_t)temp_shape_def[j]); ++ j; }
        } // for
      #endif // KERNEL2
    //#elif defined USE_MIC  // using MIC
    //  for(unsigned int i = 0; i < num_triangles * 7; ++ i)
    //    shape_def.push_back((real_t)temp_shape_def[i]);
    #else          // using CPU or MIC
//      #ifndef __SSE3__
        for(unsigned int i = 0, j = 0; i < num_triangles * CPU_T_PROP_SIZE_; ++ i) {
          if((i + 1) % CPU_T_PROP_SIZE_ == 0) shape_def.push_back((real_t) 0.0);  // padding
          else { shape_def.push_back((real_t)temp_shape_def[j]); ++ j; }
        } // for
/*      #else    // using SSE3, so store data differently: FOR CPU AND MIC (vectorization)
        #ifndef USE_MIC    // generic cpu version with SSE3 or AVX
          #ifdef INTEL_SB_AVX    // CPU version with AVX
            // group all 's', 'nx', 'ny', 'nz', 'x', 'y', 'z' together
            // for alignment at 32 bytes, make sure each of the 7 groups is padded
            // compute amount of padding
            // 32 bytes = 8 floats or 4 doubles. FIXME: assuming float only for now ...
            unsigned int padding = (8 - (num_triangles & 7)) & 7;
            unsigned int shape_size = (num_triangles + padding) * 7;
            shape_def = (real_t*) _mm_malloc(shape_size * sizeof(real_t), 32);
            if(shape_def == NULL) {
              std::cerr << "error: failed to allocate aligned memory for shape_def"
                    << std::endl;
              return 0;
            } // if
            memset(shape_def, 0, shape_size * sizeof(real_t));
            for(int i = 0; i < num_triangles; ++ i) {
              for(int j = 0; j < 7; ++ j) {
                shape_def[(num_triangles + padding) * j + i] = temp_shape_def[7 * i + j];
              } // for
            } // for
          #else        // CPU version with SSE3
            // group all 's', 'nx', 'ny', 'nz', 'x', 'y', 'z' together
            // for alignment at 16 bytes, make sure each of the 7 groups is padded
            // compute amount of padding
            // 16 bytes = 4 floats or 2 doubles. FIXME: assuming float only for now ...
            unsigned int padding = (4 - (num_triangles & 3)) & 3;
            unsigned int shape_size = (num_triangles + padding) * 7;
            shape_def = (real_t*) _mm_malloc(shape_size * sizeof(real_t), 16);
            if(shape_def == NULL) {
              std::cerr << "error: failed to allocate aligned memory for shape_def"
                    << std::endl;
              return 0;
            } // if
            memset(shape_def, 0, shape_size * sizeof(real_t));
            for(int i = 0; i < num_triangles; ++ i) {
              for(int j = 0; j < 7; ++ j) {
                shape_def[(num_triangles + padding) * j + i] = temp_shape_def[7 * i + j];
              } // for
            } // for
          #endif
        #else  // optimized for MIC only: AVX2, 64 byte alignments (512-bit vector registers)
            // FIXME: float only for now: 16 floats in one vector!
          unsigned int padding = (16 - (num_triangles & 15)) & 15;
          unsigned int shape_size = (num_triangles + padding) * 7;
          shape_def = (real_t*) _mm_malloc(shape_size * sizeof(real_t), 64);
          if(shape_def == NULL) {
            std::cerr << "error: failed to allocate aligned memory for shape_def"
                  << std::endl;
            return 0;
          } // if
          memset(shape_def, 0, shape_size * sizeof(real_t));
          for(int i = 0; i < num_triangles; ++ i) {
            for(int j = 0; j < 7; ++ j) {
              shape_def[(num_triangles + padding) * j + i] = temp_shape_def[7 * i + j];
            } // for
          } // for
          // TODO: try grouping 16 triangles together ...
          // that will give completely sequential memory access!
        #endif
      #endif // __SSE3__  */
    #endif // FF_NUM_GPU

    return num_triangles;
  } // NumericFormFactor::read_shapes_file()
} // namespace hig
