#include "hig_input.hpp"

using namespace std;
using namespace hig;

int main(int narg, char** args) {

	if(narg != 2) {
		cout << "Please specify input filename" << endl;
		exit(1);
	} // if

	if(HiGInput::instance().construct_input_config(args[1])) {
		HiGInput::instance().print_all();
		return 0;
	} // if
	return 1;
} // main()
