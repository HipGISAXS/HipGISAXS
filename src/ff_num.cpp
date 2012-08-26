/***
  *  $Id: ff_num.cpp 38 2012-08-09 23:01:20Z asarje $
  *
  *  Project:
  *
  *  File: ff_num.cpp
  *  Created: Jul 18, 2012
  *  Modified: Wed 18 Jul 2012 10:57:14 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "ff_num.hpp"

namespace hig {

	// something good may come here ...

    /**
     * Write a slice to output file
     */
    template<typename complex_t>
    void write_slice_to_file(complex_t *ff, int nqx, int nqy, int nqz, char* filename, int axis, int slice) {
        if(ff == NULL) return;

        std::cout << "** Writing slice to file ...";
        switch(axis) {
            case 0:
                break;
            default:
                std::cout << "Given axis slice writing not implemented yet!" << std::endl;
        } // switch

        if(slice >= nqx || slice < 0) {
            std::cout << "Given slice does not exist!" << std::endl;
            return;
        } // if

        std::ofstream slice_file;
        //char* outfilename = "output_ff.dat";
        char outfilename[50];
        sprintf(outfilename, "ff_p%d.dat", MPI::COMM_WORLD.Get_size());
        std::cout << " " << outfilename << " ";
        slice_file.open(outfilename);
        if(!slice_file.is_open()) {
            std::cerr << "Error opening slice file for writing." << std::endl;
            return;
        } // if

        for(int y = 0; y < nqy; ++ y) {
            for(int z = 0; z < nqz; ++ z) {
                unsigned long int index = nqx * nqy * z + nqx * y + slice;
                slice_file << ff[index].x << "," << ff[index].y << "\t";
            } // for z
            slice_file << std::endl;
        } // for y
        slice_file.close();
        std::cout << " done." << std::endl;
    } // write_slice_to_file()

} // namespace hig
