#include <iostream>
#include <time.h>

#include <analyzer/AnaInput.hpp>
#include <analyzer/Workflow.hpp>
#include <analyzer/ImageDataProc.hpp>
#include <analyzer/SimHipGISAXS.hpp>
#include <analyzer/Analyzer.hpp>
#include <analyzer/FitPOUNDERSAlgo.hpp>
#include <analyzer/SimTest.hpp>


/* The main for HipGISAXS Analyzer -- temporary ....
 */
int main(int narg, char** args) {

	if(narg > 4) {
		std::cout << "usage: hipgisaxs <input_config>" << std::endl;
		return 1;
	} // if

	hig::HipGISAXS psim(narg, args);
	if(!psim.construct_input(args[1])) {
		std::cerr << "failed to construct input containers" << std::endl;
		return 1;
	} // if
	if(!psim.fit_init()) {
		std::cerr << "error: failed to initialize hipgisaxs for fitting" << std::endl;
		return 1;
	} // if

	/////// TODO: all the following will go into the input construction ...

	//hig::ImageData img_data1("testmodel0.dat");
	hig::ImageData img_data1(psim.reference_data_path());
	hig::ImageDataProc img_data_set;
	img_data_set.push(img_data1);

	//hig::AnaInput inp;
	//inp.generate_random_vars(min,max,dim);

	/* Define objective function */
	//AbsoluteDifferenceError err;
	//hig::ObjFct pobj_fct(&err, &img_data1, &psim, psim.get_num_fit_params());
	hig::ObjFct pobj_fct(new hig::Dist(), &img_data1, &psim, psim.fit_param_init_vector().size());

	/* The algorithms to use for fitting  */
	hig::FitPOUNDERSAlgo pounders;
	pounders.set_obj_fct(&pobj_fct);
	//pounders.set_initial_vec(inp);
	pounders.set_initial_vec(psim.fit_param_init_vector());

	/* The worklow of algorithms to use for fitting  */
	hig::Workflow wf;
	wf.enq(&pounders);

	std::cout << "Setting up analyzer...\n" ;
	hig::Analyzer ana;
	//ana.set_title(title);
	ana.set_workflow(wf);
	ana.set_data(img_data_set);
	//ana.set_input(inp);
//	ana.set_input(psim.fit_param_init_vector());
	ana.set_obj_fct(&pobj_fct);

	std::cout << "Starting Analysis:\n" ;
	ana.analyze(narg,args);
	std::cout << "Done.\n";

	//test.print();

	//ana.print_output();

	return 0;
} // main()


//int main(int argc,char **argv){
//
//  /* Below are hard-coded inputs for test only. TODO: move to AnalysisBuilder class to be
//     read and generated from input file  */
//
//  /* The run title  */
//  hig::string_t title = "Test_run_0";
//
//  /******************TESTS**********************/
//  int dim=6;
//  int min= 10;
//  int max =20;
//  int n_par = 20;
//  int n_ver = 20;
//  float stpv=0, stpp=0;
//
//  hig::Data data;
//  hig::float_vec_t qy= data.create_vec(  -2, 2, stpp , n_par) ;
//  hig::float_vec_t qz= data.create_vec(  0, 2, stpv , n_ver);
//
//  /*
//  hig::float_vec_t X;
//  X.push_back(4);
//  X.push_back(4);*/
//  hig::SimTest test(dim,  qy, qz);
//  test.generate_rand_vars(min, max, dim);
//  test.compute();
//  //  test.print();
//  test.save("testmodel0.dat");
//
//  /* The forward simulation model (HipGISAXS) */
//  hig::string_t higfile = "test.anasph.hig";
//  //hig::SimHipGISAXS* psim = new hig::SimHipGISAXS(higfile);
//  hig::SimTest* psim = new hig::SimTest(dim,  qy, qz);
//
//  /* The ref. image set to fit  */
//  hig::ImageData img_data1("testmodel0.dat");
//  //  hig::ImageData img_data2("file2.dat");
//  hig::ImageDataProc img_data_set;
//  img_data_set.push(img_data1);
//  //img_data_set.push(img_data2);
//  //img_data_set.print();
//
//
//  /* The fitting parameters to fit   */
//  //hig::StructVar var1(2, "shape:s1_radius:mean");
//  //hig::StructVar var2(8, "shape:s1_radius:sigma");
//  hig::AnaInput inp;
//  inp.generate_random_vars(min,max,dim);
//  //  inp.print();
//  /******************************/
//
//  /* Define objective function */
//  hig::ObjFct* pobj_fct = new hig::ObjFct(new hig::Dist(), &img_data1, psim, dim);
//
//  /* The algorithms to use for fitting  */
//  //hig::AnaAlgorithm algo; //fit_brute_force
//
//  hig::FitPOUNDERSAlgo* pounders=new hig::FitPOUNDERSAlgo();
//  pounders->set_obj_fct(pobj_fct);
//  pounders->set_initial_vec(inp);
//
//  /*   TEST: distance computation  */
//  /*
//  hig::float_vec_t X0;
//  hig::frame_t frm;
//  X0.push_back(5.34);
//  X0.push_back(9.98);
//  pobj_fct->compute(X0);
//  hig::ImageData img_sim= pobj_fct->get_sim_data() ;
//  img_sim.print();
//  std::cout <<"diff= \n";
//  hig::float_mat_t err= pobj_fct->get_f_x();
//  hig::ImageData img_err(err, qy, qz,frm);
//  img_err.print();
//  */
//		    /*******************************/
//
//
//  /* The worklow of algorithms to use for fitting  */
//  hig::Workflow wf;
//  wf.enq(pounders);
//  //  wf.enq(pounders2);
//  //wf.enq("bruteforce");
//  // wf.print();
//
//
//  std::cout << "Setting up analyzer...\n" ;
//  hig::Analyzer ana;
//  ana.set_title(title);
//  ana.set_workflow(wf);
//  ana.set_data(img_data_set);
//  ana.set_input(inp);
//  ana.set_obj_fct(pobj_fct);
//  /*************************************/
//
//  std::cout << "Starting Analysis:\n" ;
//  ana.analyze(argc,argv);
//  std::cout << "Done.\n";
//
//  test.print();
//
//  ana.print_output();
//
//  return 0;
//}
