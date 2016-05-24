/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: input.hpp
 *
 *  Author: Dinesh Kumar
 *  Email:  dkumar@lbl.gov
 *  Date create:  
 *  Date modified: 
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

#ifndef INPUT__H
#define INPUT__H

#include <common/typedefs.hpp>
#include <model/shape.hpp>
#include <model/layer.hpp>
#include <model/unitcell.hpp>
#include <model/structure.hpp>
#include <model/inst_scattering.hpp>
#include <model/inst_detector.hpp>
#include <model/compute_params.hpp>
#include <model/fitting_params.hpp>


namespace hig {
  class Input {
    private:
      shape_list_t shapes_;
      layer_list_t layers_;
      layer_key_t layer_key_map_;
      unitcell_list_t unitcells_;
      structure_list_t structures_;
      ScatteringParams scattering_;
      DetectorParams detector_;
      ComputeParams compute_;
      FittingParams fitting_;
      
    public:
      virtual bool construct_input_config(const char *) { return false; }
      virtual bool update_params(map_t & params) { return false; }

      virtual const shape_list_t & shapes() const { return shapes_; }
      virtual const layer_list_t & layers() const { return layers_;}
      virtual const unitcell_list_t & unitcells() const { return unitcells_; }
      virtual const structure_list_t & structures() const { return structures_; }
      virtual const ScatteringParams & scattering() const { return scattering_; }
      virtual const DetectorParams & detector() const { return detector_; }
      virtual const ComputeParams & compute() const { return compute_; }
      virtual const FittingParams & fitting() const { return fitting_; }
  };
} // namespace
#endif // INPUT__H
