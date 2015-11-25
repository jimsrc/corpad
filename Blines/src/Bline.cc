#include "Bline.h"
#include "../../src/defs_turb.cc"

// DX: ancho total de la sheath
void Bline::set_Bmodel(char* fnameTURB){
	// construyo parametros del modelo, para c/region
	m.build(fnameTURB);
	// inicializar vectores I/O
	ss 	= std::vector<double>(3);
	sini 	= std::vector<double>(3);
	Lc	= m.p_turb.Lc_slab;	// [cm] long correl
}

double Bline::retLambdaMin(){
	return m.p_turb.lambda_min;	// [cm]
}

double Bline::retSigmaBoRatio(){
	return m.p_turb.sigma_Bo_ratio;	// [1]
}

void Bline::next_Brealization(){
	m.next_B_realization();
}

void Bline::set_Blineparams(double DS, int EVERY, char* fname_out, int npp){
	every	= EVERY;	// [1] nro de puntos q quiero en output, por linea de campo
	dss 	= DS*AU_in_cm;	// [cm] precision para calculo de Bline
	np 	= npp;
	//ofile.open(fname_out);
}

std::vector<std::vector<double> > Bline::returnMat(){
	return ma;
}

void Bline::calc_ds(){
	Bmod = pow(m.B[0]*m.B[0] + m.B[1]*m.B[1] + m.B[2]*m.B[2], .5); // [G] ;
        ds[0] = dss * m.B[0] / Bmod;
        ds[1] = dss * m.B[1] / Bmod;
        ds[2] = dss * m.B[2] / Bmod;
}

void Bline::init_mat(int nx, int ny){
	ma = std::vector<std::vector<double> > (nx);
	for(int i=0; i<nx; i++){
		ma[i].resize(ny);
	}
}

void Bline::calc_Bline(std::vector<double> pos){
	//double stot=0.0;
	for(i=0; i<pos.size(); i++){
		printf(" ---> pos[%d] [AU]: %g\n", i, pos[i]);
		s[i] = pos[i]*AU_in_cm;		// [cm]
	};
	std::cout << endl;
	init_mat(np, 7);

	iter 	= 0;
	nl	= 0;
	//while(abs(s[2]) < 10.*Lc){
	while(nl < np){
		//printf(" llegue auqi...\n");
		m.calc_B(s);	// calcula 'm.B[:]'
		//printf(" llegue auqi...!!!!!\n");
		calc_ds();		// valores para 'ds'

		if(iter%every==0){
			ma[nl][0] = s[0]/AU_in_cm;
			ma[nl][1] = s[1]/AU_in_cm;
			ma[nl][2] = s[2]/AU_in_cm;
			ma[nl][3] = m.B[0];
			ma[nl][4] = m.B[1];
			ma[nl][5] = m.B[2];
			ma[nl][6] = Bmod;
			nl++;
			//getchar();
		}

		for(i=0; i<3; i++){
			s[i] += ds[i];
		}
		iter++;
	}
	//stot = iter*dss;
	//printf(" stot[AU]: %g\n", stot/AU_in_cm);
	//ofile.close();
}

void Bline::report(){
	printf(" every:		%d\n", every);
	printf(" dss [AU]:	%g\n", dss/AU_in_cm);
	printf(" Lc [AU]:	%g\n", Lc/AU_in_cm);
}
