#include "read_Bmodel.h"
#include "../../src/defs_turb.cc"

void Bmodel::buildB(double DXX, char* fname_sw, char* fname_sh, char* fname_mc){
	build(DXX, fname_sw, fname_sh, fname_mc);
}

/*Bmodel::Bmodel(double DXX, char* fname_sw, char* fname_sh, char* fname_mc):
        MODEL_COMPOS(DXX, fname_sw, fname_sh, fname_mc) {}*/

void Bmodel::calc_B(std::vector<double> poss){
	double pos[3];
	pos[0] = poss[0];
	pos[1] = poss[1];
	pos[2] = poss[2];

	calc_Bfield(pos);

	BB = std::vector<double>(3);
	BB[0] = B[0];
	BB[1] = B[1];
	BB[2] = B[2];
}

std::vector<double> Bmodel::Bfield(){
	return BB;
}
