#ifndef DEFS_TURB_H
#define DEFS_TURB_H
#include "general.h"
/*-------------------------------------------------------*/

class PARAMS_SEM{
	public:
		PARAMS_SEM(void) {};
		~PARAMS_SEM(void) {};
		long slab[3];
		long two[2];
};

/*-------------------------- fases random -------------------------------*/
class FASES{
	private:
		int n_modos;

	public:
		FASES(void) {};
		//~FASES(void);
		void build(int, PARAMS_SEM *);
		double *phi_s, *a_s, *b_s;
		double *phi_2d, *b_2d;
		void construir_fases_random(PARAMS_SEM);
};

/*-------------------------- parametros turbulencia --------------------------------*/
class PARAMS_TURB{
	private:
		void read_params(string);
		void build_sigmas(void);
		void build_k_and_dk(void);
		void build_Bk_SLAB(void);
		void build_Bk_2D(void);
		void report(void);
	public:
		PARAMS_TURB(string);	// constructor
		PARAMS_TURB(void);			// constructor (otro)
		//~PARAMS_TURB(void);			// descrtructor

		string FNAME_INPUT;

		int n_modos;
		double lambda_min;
	       	double lambda_max;
		double Lc_slab, Lc_2d;	// longitudes de correlacion
		double sigma_Bo_ratio;
		double percent_slab;
		double percent_2d;

		double gS;		// potencia espectral slab
		double g2D;		// potencia espectral 2D
		double Bo;
		double sigma_S;
		double sigma_2D;
		double *dk;
		double *k;
		double *Bk_SLAB;
		double *Bk_2D;

		PARAMS_SEM sem;
		FASES fases;

		void build(string);			// puedo usarlo si es q use el contructor con "void"
};

/*-------------------------- parametros turbulencia --------------------------------*/
class MODEL_TURB{
	/*private:
		PARAMS_TURB p_turb;*/
	public:
		string FNAME_INPUT;
		MODEL_TURB(string fname_input) {build(fname_input);};			// constructor
		MODEL_TURB(void) {};			// constructor "trivial"
		//~MODEL_TURB(void);			// destructor

		void build(string);
		void calc_dB_SLAB(double *);
		void calc_dB_2D(double *);
		void calc_dB(double *);
		void calc_B(double *);

		PARAMS_TURB params_turb(void);

		double *B;		// [G]
		double *dB;		// [G]
		double *dB_SLAB;	// [G]
		double *dB_2D;		// [G]
		PARAMS_TURB p_turb;
		void next_B_realization(void);
        void fix_B_realization(int); // fija la realizacion en funcion del argumento
};

#endif //DEFS_TURB_H
