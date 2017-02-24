/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: token_mapper.hpp
 *  Created: Jun 05, 2012
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

#ifndef __TOKEN_MAPPER_HPP__
#define __TOKEN_MAPPER_HPP__

#include <unordered_map>
#include <string>

#include <config/tokens.hpp>
#include <common/enums.hpp>

namespace hig {

  class TokenMapper {

    private:

      std::unordered_map <std::string, TokenType>             KeyWords_;
      std::unordered_map <std::string, ShapeName>             ShapeKeyWords_;
      std::unordered_map <std::string, ShapeParamType>        ShapeParamKeyWords_;
      std::unordered_map <std::string, StatisticType>         StatisticKeyWords_;
      std::unordered_map <std::string, LatticeType>           LatticeKeyWords_;
      std::unordered_map <std::string, OutputRegionType>      OutputRegionKeyWords_;
      std::unordered_map <std::string, StructCorrelationType> StructCorrelationKeyWords_;
      std::unordered_map <std::string, FittingAlgorithmName>  FittingAlgorithmKeyWords_;
      std::unordered_map <std::string, FitAlgorithmParamType> FitAlgorithmParamKeyWords_;
      std::unordered_map <std::string, FittingDistanceMetric> FittingDistanceMetricKeyWords_;
      std::unordered_map <std::string, bool>                  BooleanKeyWords_;

    public:

      // this is a singleton class

      static TokenMapper& instance() {
        static TokenMapper token_mapper;
        return token_mapper;
      } // instance()

      
      TokenType get_keyword_token(const std::string& str) {
        if(KeyWords_.count(str) > 0) return KeyWords_[str];
        else return error_token;
      } // get_keyword_token()

      
      ShapeName get_shapename_token(const std::string& str) {
        if(ShapeKeyWords_.count(str) > 0) return ShapeKeyWords_[str];
        else if(has_extension(str, std::string(".shp")) ||
            has_extension(str, std::string(".hd5")) ||
            has_extension(str, std::string(".dat")) ||
            has_extension(str, std::string(".obj")))
          return shape_custom;
        else
          return shape_error;
      } // get_shapename_token()

      
      ShapeParamType get_shapeparam_token(const std::string& str) {
        if(ShapeParamKeyWords_.count(str) > 0) return ShapeParamKeyWords_[str];
        else return param_error;
      } // get_shapeparam_token()

      
      StatisticType get_stattype_token(const std::string& str) {
        if(StatisticKeyWords_.count(str) > 0) return StatisticKeyWords_[str];
        else return stat_error;
      } // get_stattype_token()


      LatticeType get_lattice_type(const std::string& str) {
        if(LatticeKeyWords_.count(str) > 0) return LatticeKeyWords_[str];
        else return lattice_error;
      } // get_lattice_type_token()


      OutputRegionType get_output_region_type(const std::string& str) {
        if(OutputRegionKeyWords_.count(str) > 0) return OutputRegionKeyWords_[str];
        else return region_error;
      } // get_output_region_type()


      StructCorrelationType get_structcorr_type(const std::string& str) {
        if(StructCorrelationKeyWords_.count(str) > 0) return StructCorrelationKeyWords_[str];
        else return structcorr_error;
      } // get_structcorr_type()

      
      bool key_exists(const std::string& str) {
        if(KeyWords_.count(str) > 0 ||
            ShapeKeyWords_.count(str) > 0 ||
            ShapeParamKeyWords_.count(str) > 0 ||
            StatisticKeyWords_.count(str) > 0 ||
            LatticeKeyWords_.count(str) > 0 ||
            OutputRegionKeyWords_.count(str) > 0)
          return true;
        return false;
      } // keyword_exists()


      bool keyword_token_exists(TokenType token) {
        // a sequential search for value 'token'
        std::unordered_map <std::string, TokenType>::iterator i = KeyWords_.begin();
        while(i != KeyWords_.end()) { if((*i).second == token) return true; i ++; }
        return false;
      } // token_exists()


      OutputRegionType get_reference_data_region_type(std::string str) {
        if(OutputRegionKeyWords_.count(str) > 0) return OutputRegionKeyWords_[str];
        else return region_error;
      } // get_reference_data_region_type()


      FitAlgorithmParamType get_fit_algorithm_param_token(std::string str) {
        if(FitAlgorithmParamKeyWords_.count(str) > 0) return FitAlgorithmParamKeyWords_[str];
        else return algo_param_error;
      } // get_fit_algorithm_param_token()


      FittingAlgorithmName get_fit_algorithm_name(std::string str) {
        if(FittingAlgorithmKeyWords_.count(str) > 0) return FittingAlgorithmKeyWords_[str];
        else return algo_error;
      } // get_fit_algorithm_name()


      FittingDistanceMetric get_fit_distance_metric(std::string str) {
        if(FittingDistanceMetricKeyWords_.count(str) > 0) return FittingDistanceMetricKeyWords_[str];
        else return metric_error;
      } // 


      bool get_boolean(std::string str) {
        if(BooleanKeyWords_.count(str) > 0) return BooleanKeyWords_[str];
        else return false;
      } // get_boolean()


    private:

      /* constructor */
      TokenMapper() {

        /* language keywords */

        KeyWords_[std::string("a")]               = struct_grain_lattice_a_token;
        KeyWords_[std::string("abangle")]         = struct_grain_lattice_abangle_token;
        KeyWords_[std::string("algorithm")]       = fit_algorithm_token;
        KeyWords_[std::string("algoname")]        = fit_algorithm_name_token;
        KeyWords_[std::string("algoorder")]       = fit_algorithm_order_token;
        KeyWords_[std::string("algoparam")]       = fit_algorithm_param_token;
        KeyWords_[std::string("alphai")]          = instrument_scatter_alphai_token;
        KeyWords_[std::string("angles")]          = struct_ensemble_orient_rot_angles_token;
        KeyWords_[std::string("anglelocation")]   = struct_ensemble_orient_rot_anglelocation_token;
        KeyWords_[std::string("anglemean")]       = struct_ensemble_orient_rot_anglemean_token;
        KeyWords_[std::string("anglescale")]      = struct_ensemble_orient_rot_anglescale_token;
        KeyWords_[std::string("anglesd")]         = struct_ensemble_orient_rot_anglesd_token;
        KeyWords_[std::string("axis")]            = struct_ensemble_orient_rot_axis_token;
        KeyWords_[std::string("b")]               = struct_grain_lattice_b_token;
        KeyWords_[std::string("beta")]            = refindex_beta_token;
        KeyWords_[std::string("c")]               = struct_grain_lattice_c_token;
        KeyWords_[std::string("caratio")]         = struct_grain_lattice_caratio_token;
        KeyWords_[std::string("coherence")]       = instrument_scatter_coherence_token;
        KeyWords_[std::string("computation")]     = compute_token;
        KeyWords_[std::string("delta")]           = refindex_delta_token;
        KeyWords_[std::string("detector")]        = instrument_detector_token;
        KeyWords_[std::string("dimensions")]      = struct_dims;
        KeyWords_[std::string("directbeam")]      = instrument_detector_dirbeam_token;
        KeyWords_[std::string("distancemetric")]  = fit_algorithm_distance_metric_token;
        KeyWords_[std::string("distribution")]    = struct_ensemble_distribution_token;
        KeyWords_[std::string("domain")]          = struct_paracrystal_domain_size;
        KeyWords_[std::string("element")]         = unitcell_element_token;
        KeyWords_[std::string("ensemble")]        = struct_ensemble_token;
        KeyWords_[std::string("expt")]            = instrument_scatter_expt_token;
        KeyWords_[std::string("fitparam")]        = fit_param_token;
        KeyWords_[std::string("fitregion")]       = fit_reference_data_region_token;
        KeyWords_[std::string("fitting")]         = fit_token;
        KeyWords_[std::string("grain")]           = struct_grain_token;
        KeyWords_[std::string("hipGisaxsInput")]  = hipgisaxs_token;
        KeyWords_[std::string("hkl")]             = struct_grain_lattice_hkl_token;
        KeyWords_[std::string("include")]         = include_token;
        KeyWords_[std::string("init")]            = fit_param_init_token;
        KeyWords_[std::string("inplanerot")]      = instrument_scatter_inplanerot_token;
        KeyWords_[std::string("instrumentation")] = instrument_token;
        KeyWords_[std::string("iratio")]          = struct_iratio_token;
        KeyWords_[std::string("key")]             = key_token;
        KeyWords_[std::string("lattice")]         = struct_grain_lattice_token;
        KeyWords_[std::string("layer")]           = layer_token;
        KeyWords_[std::string("layer:key")]       = struct_grain_lkey_token;
        KeyWords_[std::string("locations")]       = unitcell_element_locations_token;
        KeyWords_[std::string("mask")]            = fit_reference_data_mask_token;
        KeyWords_[std::string("max")]             = max_token;
        KeyWords_[std::string("maxgrains")]       = struct_ensemble_maxgrains_token;
        KeyWords_[std::string("maxpoint")]        = compute_outregion_maxpoint_token;
        KeyWords_[std::string("mean")]            = mean_token;
        KeyWords_[std::string("method")]          = compute_method_token;
        KeyWords_[std::string("min")]             = min_token;
        KeyWords_[std::string("minpoint")]        = compute_outregion_minpoint_token;
        KeyWords_[std::string("name")]            = shape_name_token;
        KeyWords_[std::string("npoints")]         = fit_reference_data_npoints_token;
        KeyWords_[std::string("nsamples")]        = nsamples_token;
        KeyWords_[std::string("nslices")]         = compute_nslices_token;
        KeyWords_[std::string("nvalues")]         = shape_param_nvalues_token;
        KeyWords_[std::string("order")]           = layer_order_token;
        KeyWords_[std::string("orientations")]    = struct_ensemble_orient_token;
        KeyWords_[std::string("origin")]          = instrument_detector_origin_token;
        KeyWords_[std::string("originvec")]       = shape_originvec_token;
        KeyWords_[std::string("outputregion")]    = compute_outregion_token;
        KeyWords_[std::string("p1")]              = shape_param_p1_token;    // mean
        KeyWords_[std::string("p2")]              = shape_param_p2_token;    // std dev
        KeyWords_[std::string("palette")]         = compute_palette_token;
        KeyWords_[std::string("paracrystal")]     = struct_paracrystal;
        KeyWords_[std::string("parallel")]        = fit_reference_data_npoints_parallel_token;
        KeyWords_[std::string("param")]           = shape_param_token;
        KeyWords_[std::string("path")]            = fit_reference_data_path_token;
        KeyWords_[std::string("pathprefix")]      = compute_path_token;
        KeyWords_[std::string("percusyevick")]    = struct_percusyevick;
        KeyWords_[std::string("perpendicular")]   = fit_reference_data_npoints_perpendicular_token;
        KeyWords_[std::string("photon")]          = instrument_scatter_photon_token;
        KeyWords_[std::string("pixelsize")]       = instrument_detector_pixsize_token;
        KeyWords_[std::string("polarization")]    = instrument_scatter_polarize_token;
        KeyWords_[std::string("pvalue")]          = fit_algorithm_param_value_token;
        KeyWords_[std::string("range")]           = fit_param_range_token;
        KeyWords_[std::string("regmax")]          = fit_reference_data_region_max_token;
        KeyWords_[std::string("regmin")]          = fit_reference_data_region_min_token;
        KeyWords_[std::string("regularization")]  = fit_algorithm_regularization_token;
        KeyWords_[std::string("repetition")]      = struct_grain_repetition_token;
        KeyWords_[std::string("repetitiondist")]  = struct_grain_repetitiondist_token;
        KeyWords_[std::string("restart")]         = fit_algorithm_restart_token;
        KeyWords_[std::string("referencedata")]   = fit_reference_data_token;
        KeyWords_[std::string("refindex")]        = refindex_token;
        KeyWords_[std::string("resolution")]      = compute_resolution_token;
        KeyWords_[std::string("rot")]             = rot_token;
        KeyWords_[std::string("rot1")]            = struct_ensemble_orient_rot1_token;
        KeyWords_[std::string("rot2")]            = struct_ensemble_orient_rot2_token;
        KeyWords_[std::string("rot3")]            = struct_ensemble_orient_rot3_token;
        KeyWords_[std::string("runname")]         = compute_runname_token;
        KeyWords_[std::string("saveff")]          = compute_saveff_token;
        KeyWords_[std::string("savesf")]          = compute_savesf_token;
        KeyWords_[std::string("scaling")]         = struct_grain_scaling_token;
        KeyWords_[std::string("scattering")]      = instrument_scatter_token;
        KeyWords_[std::string("sdd")]             = instrument_detector_sdd_token;
        KeyWords_[std::string("shape")]           = shape_token;
        KeyWords_[std::string("shape:key")]       = unitcell_element_skey_token;
        KeyWords_[std::string("smearing")]        = instrument_scatter_smearing_token;
        KeyWords_[std::string("spacing")]         = struct_ensemble_spacing_token;
        KeyWords_[std::string("spotarea")]        = instrument_scatter_spotarea_token;
        KeyWords_[std::string("stat")]            = stat_token;
        KeyWords_[std::string("stddev")]          = stddev_token;
        KeyWords_[std::string("step")]            = step_token;
        KeyWords_[std::string("structcorrelation")] = compute_structcorr_token;
        KeyWords_[std::string("structure")]       = struct_token;
        KeyWords_[std::string("thickness")]       = layer_thickness_token;
        KeyWords_[std::string("tilt")]            = instrument_scatter_tilt_token;
        KeyWords_[std::string("tolerance")]       = fit_algorithm_tolerance_token;
        KeyWords_[std::string("totalpixels")]     = instrument_detector_totpix_token;
        KeyWords_[std::string("transvec")]        = struct_grain_transvec_token;
        KeyWords_[std::string("type")]            = type_token;
        KeyWords_[std::string("unit")]            = instrument_scatter_photon_unit_token;
        KeyWords_[std::string("unitcell")]        = unitcell_token;
        KeyWords_[std::string("unitcell:key")]    = struct_grain_ukey_token;
        KeyWords_[std::string("value")]           = instrument_scatter_photon_value_token;
        KeyWords_[std::string("variable")]        = fit_param_variable_token;
        KeyWords_[std::string("volfraction")]     = struct_percusyevick_volfract;
        KeyWords_[std::string("xrepetition")]     = struct_grain_xrepetition_token;
        KeyWords_[std::string("xrot")]            = shape_xrot_token;
        KeyWords_[std::string("xspacing")]        = struct_paracrystal_xspacing;
        KeyWords_[std::string("yrepetition")]     = struct_grain_yrepetition_token;
        KeyWords_[std::string("yrot")]            = shape_yrot_token;
        KeyWords_[std::string("yspacing")]        = struct_paracrystal_yspacing;
        KeyWords_[std::string("zrepetition")]     = struct_grain_zrepetition_token;
        KeyWords_[std::string("zrot")]            = shape_zrot_token;
      
        /* shape name keywords */

        ShapeKeyWords_[std::string("box")]            = shape_box;
        ShapeKeyWords_[std::string("cube")]           = shape_cube;
        ShapeKeyWords_[std::string("cylinder")]       = shape_cylinder;
        ShapeKeyWords_[std::string("hcylinder")]      = shape_horizontal_cylinder;
        ShapeKeyWords_[std::string("randcylinders")]  = shape_random_cylinders;
        ShapeKeyWords_[std::string("sphere")]         = shape_sphere;
        ShapeKeyWords_[std::string("pyramid")]        = shape_pyramid;
        ShapeKeyWords_[std::string("trunccone")]      = shape_trunccone;
        ShapeKeyWords_[std::string("prism3")]         = shape_prism3;
        ShapeKeyWords_[std::string("prism6")]         = shape_prism6;
        ShapeKeyWords_[std::string("prism3x")]        = shape_prism3x;
        ShapeKeyWords_[std::string("sawtooth")]       = shape_sawtooth;
        ShapeKeyWords_[std::string("custom")]         = shape_custom;
      
        /* shape parameter type keywords */

        ShapeParamKeyWords_[std::string("radius")]    = param_radius;
        ShapeParamKeyWords_[std::string("xsize")]     = param_xsize;
        ShapeParamKeyWords_[std::string("ysize")]     = param_ysize;
        ShapeParamKeyWords_[std::string("height")]    = param_height;
        ShapeParamKeyWords_[std::string("edge")]      = param_edge;
        ShapeParamKeyWords_[std::string("baseangle")] = param_baseangle;
      
        /* statistic type keywords */

        StatisticKeyWords_[std::string("none")]     = stat_none;
        StatisticKeyWords_[std::string("single")]   = stat_none;
        StatisticKeyWords_[std::string("uniform")]  = stat_uniform;
        StatisticKeyWords_[std::string("random")]   = stat_random;
        StatisticKeyWords_[std::string("range")]    = stat_range;
        StatisticKeyWords_[std::string("gaussian")] = stat_gaussian;
        StatisticKeyWords_[std::string("cauchy")]   = stat_cauchy;
        StatisticKeyWords_[std::string("t")]        = stat_t;

        /* lattice type keywords */

        LatticeKeyWords_[std::string("bcc")]    = lattice_bcc;
        LatticeKeyWords_[std::string("cubic")]  = lattice_cubic;
        LatticeKeyWords_[std::string("fcc")]    = lattice_fcc;
        LatticeKeyWords_[std::string("fco")]    = lattice_fco;
        LatticeKeyWords_[std::string("hcp")]    = lattice_hcp;
        LatticeKeyWords_[std::string("hex")]    = lattice_hex;

        /* output region type keywords */

        OutputRegionKeyWords_[std::string("angles")]  = region_angles;
        OutputRegionKeyWords_[std::string("pixels")]  = region_pixels;
        OutputRegionKeyWords_[std::string("qspace")]  = region_qspace;

        /* structure grain/ensemble correlation keywords */

        StructCorrelationKeyWords_[std::string("nGnE")] = structcorr_nGnE;
        StructCorrelationKeyWords_[std::string("nGE")]  = structcorr_nGE;
        StructCorrelationKeyWords_[std::string("GnE")]  = structcorr_GnE;
        StructCorrelationKeyWords_[std::string("GE")]   = structcorr_GE;

        /* fitting algorithm name keywords */

        FittingAlgorithmKeyWords_[std::string("lmvm")]          = algo_lmvm;
        FittingAlgorithmKeyWords_[std::string("pounders")]      = algo_pounders;
        FittingAlgorithmKeyWords_[std::string("pso")]           = algo_pso;
        FittingAlgorithmKeyWords_[std::string("bruteforce")]    = algo_bruteforce;
        FittingAlgorithmKeyWords_[std::string("none_pounders")] = algo_none_pounders;

        /* fitting algorithm parameter keywords */

        // pounders
        FitAlgorithmParamKeyWords_[std::string("pounders_delta")]       = algo_pounders_param_delta;

        // pso
        FitAlgorithmParamKeyWords_[std::string("pso_omega")]            = algo_pso_param_omega;
        FitAlgorithmParamKeyWords_[std::string("pso_phi1")]             = algo_pso_param_phi1;
        FitAlgorithmParamKeyWords_[std::string("pso_phi2")]             = algo_pso_param_phi2;
        FitAlgorithmParamKeyWords_[std::string("pso_num_particles")]    = algo_pso_param_nparticle;
        FitAlgorithmParamKeyWords_[std::string("pso_num_generations")]  = algo_pso_param_ngen;
        FitAlgorithmParamKeyWords_[std::string("pso_tune_omega")]       = algo_pso_param_tune_omega;
        FitAlgorithmParamKeyWords_[std::string("pso_type")]             = algo_pso_param_type;

        /* fitting distance metric keywords */

        FittingDistanceMetricKeyWords_[std::string("sqrt_unit_l1")]           = metric_sqrt_unit_norm_l1;
        FittingDistanceMetricKeyWords_[std::string("sqrt_unit_l2")]           = metric_sqrt_unit_norm_l2;
        FittingDistanceMetricKeyWords_[std::string("sqrt_c_l2")]              = metric_sqrt_c_norm_l2;
        FittingDistanceMetricKeyWords_[std::string("cbrt_unit_l1")]           = metric_cbrt_unit_norm_l1;
        FittingDistanceMetricKeyWords_[std::string("cbrt_unit_l2")]           = metric_cbrt_unit_norm_l2;
        FittingDistanceMetricKeyWords_[std::string("cbrt_c_l2")]              = metric_cbrt_c_norm_l2;
        FittingDistanceMetricKeyWords_[std::string("sqrt_unit_l1_residual")]  = metric_sqrt_unit_norm_l1_residual;
        FittingDistanceMetricKeyWords_[std::string("sqrt_unit_l2_residual")]  = metric_sqrt_unit_norm_l2_residual;
        FittingDistanceMetricKeyWords_[std::string("sqrt_c_l2_residual")]     = metric_sqrt_c_norm_l2_residual;

        /* boolean keywords */

        BooleanKeyWords_[std::string("true")]     = true;
        BooleanKeyWords_[std::string("True")]     = true;
        BooleanKeyWords_[std::string("TRUE")]     = true;
        BooleanKeyWords_[std::string("yes")]      = true;
        BooleanKeyWords_[std::string("Yes")]      = true;
        BooleanKeyWords_[std::string("YES")]      = true;
        BooleanKeyWords_[std::string("false")]    = false;
        BooleanKeyWords_[std::string("False")]    = false;
        BooleanKeyWords_[std::string("FALSE")]    = false;
        BooleanKeyWords_[std::string("no")]       = false;
        BooleanKeyWords_[std::string("No")]       = false;
        BooleanKeyWords_[std::string("NO")]       = false;

    } // TokenMapper()

    // singleton
    TokenMapper(const TokenMapper&);
    TokenMapper& operator=(const TokenMapper&);

    bool has_extension(const std::string& s, const std::string& e) const {
      unsigned int p = s.rfind(e);
      if(p != s.length() - e.length())  // s is not a suffix of
        return false;
      return true;
    } // has_extension()

  }; // class TokenMapper

} // namespace hig

#endif /* __TOKEN_MAPPER_HPP__ */
