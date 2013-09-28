/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: object2hdf5.h
 *  Created: Aug 25, 2012
 *  Modified: Fri 27 Sep 2013 09:02:30 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
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

#ifndef _OBJECT2HDF5_H_
#define _OBJECT2HDF5_H_

#include <hdf5.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

void h5_shape_reader(const char* hdf5_filename, double** shape_def, unsigned int* num_triangles);

#ifdef __cplusplus
}
#endif

#endif // _OBJECT2HDF5_H_
