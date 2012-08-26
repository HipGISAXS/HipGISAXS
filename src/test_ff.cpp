/***
  *  $Id: test_ff.cpp 38 2012-08-09 23:01:20Z asarje $
  *
  *  Project:
  *
  *  File: test_ff.cpp
  *  Created: Aug 01, 2012
  *  Modified: Thu 02 Aug 2012 03:44:19 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <mpi.h>
#include <cuComplex.h>

#include "qgrid.hpp"
#include "ff_num.hpp"

using namespace std;
using namespace hig;

int main(int narg, char** args) {
	MPI::Init(narg, args);

	const char* filename = args[1];
	QGrid::instance().create_test();
	NumericFormFactor<float, cuFloatComplex> nff(4);
	nff.init();
	cuFloatComplex *ff = NULL;
	nff.compute(filename, ff, MPI::COMM_WORLD);

	std::cout << "FF: " /*<< QGrid::instance().nqx() << ", "
						<< QGrid::instance().nqy() << ", "
						<< QGrid::instance().nqz_extended()*/ << std::endl;
	for(int k = 0; k < QGrid::instance().nqz_extended(); ++ k) {
		//std::cout << " + " << k << ";" << std::endl;
		for(int j = 0; j < QGrid::instance().nqy(); ++ j) {
			for(int i = 0; i < QGrid::instance().nqx(); ++ i) {
				int index = k * QGrid::instance().nqx() * QGrid::instance().nqy() +
							j * QGrid::instance().nqx() + i;
				std::cout << /*QGrid::instance().qz(k % QGrid::instance().nqz())
							<< " -> " << QGrid::instance().qz_extended(k)
							<< " -> " <<*/ ff[index].x << "," << ff[index].y << "    ";
			} // for
			std::cout << std::endl;
		} // for
		std::cout << std::endl;
	} // for
	std::cout << std::endl;

	if(ff != NULL) delete[] ff;
	MPI::Finalize();
	return 0;
} // main()
