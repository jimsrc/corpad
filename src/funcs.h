#ifndef FUNCS_H
#define FUNCS_H

//---------------------------------------------------
class PARAMS : public MODEL_TURB{
	public:
		PARAMS(string);
		void calc_Bfield(VecDoub_I &);
		//PARAMS & operator=(const PARAMS &rhs);
	private:
		double pos[3];
};

//---------------------------------------------------
template <class Stepper>
class Output {
	public:
		Int kmax;
		Int nvar;
		Int nsave;
		bool dense;
		Int count;
		Doub x1,x2,xout,dxout;
		VecDoub xsave;
		MatDoub ysave;
		void build(string, Int, Doub, Int, char*);
		Output(void);
		Output(string, const Int, char*);
		void init(const Int, const Doub, const Doub);
		void resize(void);
		void save_dense(Stepper &, const Doub, const Doub);
		void save(const Doub, VecDoub_I &);
		void out(const Int,const Doub,VecDoub_I &,Stepper &,const Doub);
		void save2file(void);
		ofstream ofile_trj, ofile_misc;
		bool file_exist(void);
		void resizeTau(void);
		//esto lo agrego para guardar cosas de la historia de 
		//las trayectorias:
		int nfilTau, ncolTau;		// tamanio para 'Tau'
		int nreb;			// nro de rebotes/scatterings en pitch
		MatDoub Tau;			// tiempo de camino libre medio paralelo, y su posic x
		VecDoub mu;
		void set_Bmodel(PARAMS);	// para apuntar al modelo q uso en main()
		void tic(void), toc(void);	// cronometro para c/pla
		Doub trun;			// tiempo de simulacion de c/pla
		Int nsteps;			// nro total de pasos de c/pla
	private:
		char fname_trj[200], fname_misc[200];	// nombres de archivos de salida
		PARAMS *pm;
		double pos[3], vmod, bmod;
		void save_pitch(void);
		double bx, by, bz, vx, vy, vz;
		VecDoub XSaveGen;	// tiempos de la salida
		int n_tscales, cc, ndec, inid, nd, base, maxd;	
		string str_tscale;	// tipo de escala temporal para la salida
		double decade, dt;
		void set_savetimes(Doub);
		//----- histo del 'Tau'
		MatDoub HistTau;
		double dTau, maxTau, avrTau;
		int nHistTau, nTau, dimHistTau;
		void build_HistTau(void);
};

#endif //FUNCS_H
