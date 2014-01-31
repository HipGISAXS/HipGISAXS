/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: ObjFct.cpp
 *  Created: Dec 26, 2013
 *  Modified: Fri 31 Jan 2014 01:57:15 PM PST
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
*
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <map>
#include <analyzer/ObjFct.hpp>

namespace hig{


  PetscErrorCode EvaluateFunction(TaoSolver tao, Vec X, Vec F, void *ptr)
  {
    // Compute F(X)
    PetscFunctionBegin;
    VecView(X, PETSC_VIEWER_STDOUT_WORLD);

    PetscErrorCode ierr;
    PetscReal *x, *f;

    ierr = VecGetArray(X,&x);// CHKERRQ(ierr);
    ierr = VecGetArray(F,&f); //CHKERRQ(ierr);

    PetscReal* ff =( (ObjFct*) ptr)->tao_compute(x);
    int nobs = ( (ObjFct*) ptr)->get_nobs();

    for(int i=0; i< nobs; i++){
      f[i]= ff[i];
    }

    ierr = VecRestoreArray(X,&x); CHKERRQ(ierr);
    ierr = VecRestoreArray(F,&ff); CHKERRQ(ierr);

    std::cout << "Eval X=\n" ;
    VecView(X, PETSC_VIEWER_STDOUT_WORLD);
    std::cout << " = \n" ;
    //VecView(F, PETSC_VIEWER_STDOUT_WORLD);

    PetscFunctionReturn(0);
    return 0;
  }

  PetscErrorCode EvaluateJacobian(TaoSolver tao, Vec X, Mat *J, Mat *Jpre, MatStructure*matstruct,void *ptr)
  {
    Mat JJ =  (( (ObjFct *)ptr )->tao_compute_jacobian(X) );
    J = &JJ;
    PetscFunctionReturn(0);
    return 0;
  }

  /*****************************************************************************/


  float_mat_t ObjFct::compute(){

    //    return {0};
  }


  float_mat_t ObjFct::compute(float_vec_t X){
    if(psim_)
      {
	//if (is_valid_ && pdata_ref_){
//	psim_->update_vars(X);
//	psim_->compute();
	//	pdata_sim_=psim_->get_data();
//	f_x_ = pdist_->dist( psim_->get_data()   , *pdata_ref_ );
	return f_x_;
      }
  }

	////// temporary .. for hipgisaxs -- abhinav
	float_mat_t ObjFct::operator()(float_vec_t x) {
		float_t *gisaxs_data = NULL;
		// construct param_vals ...
		std::vector <std::string> params = psim_->get_fit_param_keys();
		std::map <std::string, float_t> param_vals;
		for(int i = 0; i < x.size(); ++ i) {
			param_vals[params[i]] = x[i];
		} // for
		psim_->update_params(param_vals);
		psim_->compute_gisaxs(gisaxs_data);
		float_t* ref_data = (*pdata_ref_).data();
		f_x_ = (*pdist_).dist(gisaxs_data, ref_data, (*pdata_ref_).size());
		//std::cout << "@@ ERROR: " << err << std::endl;
		return f_x_;
	} // ObjFct::operator()()


  float_mat_t ObjFct::compute_jacobian(float_vec_t X, int dim){
    if(psim_)
      {
        //if (is_valid_ && pdata_ref_){
	for (int n=0; n<dim; n++){


	  //f_xpdx_n = compute(X+DXn);


	}

	/*	psim_->update_vars(X);
        psim_->compute();
        pdata_sim_=psim_->get_data();
        J_x_ = pdist_->dist(*pdata_sim_ , *pdata_ref_ );
	*/
	return J_x_;
      }
  }

  float_mat_t ObjFct::compute_test(float_vec_t X){
    if (X.size() == 2){
      f_x_.clear();
      float x0= X[0];
      float x1= X[1];
      float fx= exp(-( (x0- 0 )*(x0 - 0) + (x1 - 0)*(x1 - 0) )/ 2 );
      //float_vec_t v;
      f_x_.push_back(fx);
      //f_x_.push_back(v);
      return f_x_;
    }
  }

  PetscReal* ObjFct::tao_compute(PetscReal* x){
    /* TODO: lots of overhead!!!! improve performance!!!  */
    //PetscFunctionBegin;
    PetscReal* pfx = new PetscReal[get_nobs()];
    float_vec_t X;

    for(int i=0; i<dim_; i++){
      X.push_back(x[i]);
    }

    /* Compute F(x) */
    //compute(X);
    (*this)(X);
    /************** */

    int ix=0;
    for(int iv=0; iv<n_ver_; iv++)
      for(int ip=0; ip<n_par_; ip++){
	pfx[ix] = f_x_[ix];
	ix++;
      }
//	pfx[0] = f_x_[0][0];
    return pfx;
  }

  Mat ObjFct::tao_compute_jacobian(Vec X){

    /* TODO: lots of overhead!!!! improve performance!!!  */
    //PetscFunctionBegin;
    Vec fx;
    float_vec_t x;
    //    PetscErrorCode ierr;
    PetscInt dim;
    PetscScalar y[1]={0};

    VecGetSize(X, &dim);

    for(int i=0; i<dim; i++){
      VecGetValues(X, 1, &i , y);
      x.push_back(y[0]);
    }

    VecView(X, PETSC_VIEWER_STDOUT_WORLD);

    /* Compute J(x) */
    compute_jacobian(x, dim);
    /*   */
    /*
    ierr = VecGetArray(DX,&dx);

    ierr = VecCreateSeq(PETSC_COMM_SELF, nn, &DXI);

    ierr = VecCreateSeq(PETSC_COMM_SELF, nn, &XP1);
    ierr = VecCreateSeq(PETSC_COMM_SELF, nn, &XM1);
    // Compute G(X)
    for (i=0; i<nn; i++){

      ierr = VecSet(XP1, 0);
      ierr = VecGetArray(XP1,&xp1);

      ierr = VecSet(XM1, 0);
      ierr = VecGetArray(XM1,&xm1);

      ierr = VecSet(DXI, 0);
      VecSetValue(DXI, i, dx[i]/2, INSERT_VALUES);
      //VecView(DXI,PETSC_VIEWER_STDOUT_WORLD);

      ierr = VecGetArray(DXI, &dxi);

      for (j=0;j<nn;j++){
        xp1[j] = x[j] + dxi[j];
        xm1[j] = x[j] - dxi[j];
      }

      ierr = VecRestoreArray(XP1,&xp1);
      ierr = VecRestoreArray(XM1,&xm1);
      ierr = VecRestoreArray(DXI,&dxi);

      PetscReal* fxp1 = TestFormFunction(XP1, user);
      PetscReal* fxm1 = TestFormFunction(XM1, user);

      for (k=0; k<nobsy*nobsz; k++) {
        user->j[k][i] = (fxp1[k] - fxm1[k])/dx[i];
      }
    }
    // Assemble the matrix
    ierr = MatSetValues(*J,nobsy*nobsz,user->idm, NPARAMETERS, user->idn,(PetscReal *)user->j,
                        INSERT_VALUES); //CHKERRQ(ierr);
    ierr = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);

    //PetscPrintf(PETSC_COMM_SELF,"\n---- J(X) (Jacobian) =  -----\n");
    //MatView(*J,PETSC_VIEWER_STDOUT_WORLD);

    // Return the handles
    ierr = VecRestoreArray(X,&x); //CHKERRQ(ierr);
    //PetscLogFlops(user->n * 13);

    */
  }

}
