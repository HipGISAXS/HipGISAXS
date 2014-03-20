/**
 *	Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *	File: ImageData.cpp
 *	Created: Dec 26, 2013
 *	Modified: Thu 20 Mar 2014 02:07:35 PM PDT
 *
 *	Author: Slim Chourou <stchourou@lbl.gov>
 *	Developers: Slim Chourou <stchourou@lbl.gov>
 *							Abhinav Sarje <asarje@lbl.gov>
 *							Alexander Hexemer <ahexemer@lbl.gov>
 *							Xiaoye Li <xsli@lbl.gov>
 *
 *	Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *	used by employees of academic research institutions, not-for-profit
 *	research laboratories, or governmental research facilities. Please read the
 *	accompanying LICENSE file before downloading the software. By downloading
 *	the software, you are agreeing to be bound by the terms of this
 *	NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <analyzer/ImageData.hpp>

namespace hig {

	//float ImageData::img_p(int i, int j) const {
	float_t ImageData::operator()(int i, int j) const {
		if(i < 0 || i >= n_par_ || j < 0 || j >= n_ver_) return 0;
		return data_[j * n_par_ + i];
	} // ImageData::operator()()


	float_t ImageData::operator()(unsigned int i, unsigned int j) const {
		if(i >= n_par_ || j >= n_ver_) return 0;
		return data_[j * n_par_ + i];
	} // ImageData::operator()()


	float_t ImageData::operator()(float_t i, float_t j) const {
		// TODO ...
		return 0;
	} // ImageData::operator()()

	//float ImageData::img_q(float qv, float qp) const {
	//	return 0;
	//}

	void	ImageData::print() const {
		std::cout << filename_ << std::endl;
		for(int iv=0; iv<n_ver_; iv++)
			{
	for(int ip=0; ip<n_par_; ip++)
		{
			std::cout << data_[iv * n_par_ + ip] << "	";
		}
	std::cout <<std::endl;
			}
	}

	void	ImageData::save(string_t filename) const {
		std::ofstream file;
		file.open(filename.c_str());

		for(int iv=0; iv<n_ver_; iv++)
			{
				for(int ip=0; ip<n_par_; ip++)
					{
			file << data_[iv * n_par_ + ip] << "	";
					}
	file << "\n";
			}
		file.close();
	}

	float_vec_t ImageData::read_string_values(string_t line){
		float_vec_t array;

		std::stringstream ssin(line);
		std::copy(	std::istream_iterator<float>(ssin),
		std::istream_iterator<float>(),
		std::back_inserter(array));

		return array;
	}

	void ImageData::read(string_t filename){
		int nv = -1;
		int np = 0;
		data_.clear();
		std::string line;
		std::ifstream file(filename.c_str());
		if (file.is_open()) { //if the file is open
			while (!file.eof()) { //while the end of file is NOT reached
				getline(file,line); //get one line from the file
				nv++;
				float_vec_t img_z = read_string_values(line);
				if(nv == 0) np = img_z.size();
				data_.insert(data_.end(), img_z.begin(), img_z.end());
			} // while
			file.close(); //closing the file
			n_par_ = np;
			n_ver_ = nv;
			//convert_data();
		} // if
		else std::cerr << "error: unable to open file " << filename << std::endl;
	} // ImageData::read()


	// temporary .... -- abhinav
	//float_t* ImageData::convert_data() {
	//	std::cout << "ALLOCATING DATA MEMORY" << std::endl;
	//	data_.clear();
	//	unsigned int i = 0;
	//	for(float_mat_t::iterator r = img_.begin(); r != img_.end(); ++ r) {
	//		//for(float_vec_t::iterator c = (*r).begin(); c != (*r).end(); ++ c) {
	//			data_.push_back(*r);
	//		//} // for
	//	} // for
	//} // ImageData::data()

}
