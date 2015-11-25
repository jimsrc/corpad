//#include "../../src/general.cc"
#include "../../src/defs_turb.h"

class Bmodel : public MODEL_COMPOS{
        public:
		Bmodel(void){};
                //Bmodel(double, char*, char*, char*);
		void buildB(double, char*, char*, char*);
		void calc_B(std::vector<double>);
		std::vector<double> BB;
		std::vector<double> Bfield();
};
