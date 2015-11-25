#include "../../src/defs_turb.h"

class Bline{
	public:
		Bline(void) {};
		//~Bline(void) {};
		//--- en vez de una sola funcion "build", hago:
		void set_Bmodel(char*);	// 
		void set_Blineparams(double,int,char*,int);		// 
		//--------------------------------------------
		std::vector<double> ss, sini;
		void calc_Bline(std::vector<double>); // recbie posicion en [cm]
		int every;
		void report(void);
		std::vector<std::vector<double> > returnMat(void);
		double retLambdaMin(void);
		void next_Brealization(void);
		double retSigmaBoRatio(void);
	private:
		MODEL_TURB m;
		double s[3], ds[3], dss;
		double Lc;		// [cm] long correl del SW
		double Bmod;		// [G]
		int iter;		// nro de iteraciones
		int i, nl, np;
		void calc_ds(void);
		ofstream ofile;
		std::vector<std::vector<double> > ma;
		void init_mat(int, int);
};
