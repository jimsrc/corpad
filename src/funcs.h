#ifndef FUNCS_H
#define FUNCS_H

#ifdef CYTHON
    #define PRIVATE_OR_PUBLIC public
#else
    #define PRIVATE_OR_PUBLIC private
#endif //CYTHON

//#include "control.h"
//#include "general.h"
#include "nr3.h"
#include "defs_turb.h"

//---------------------------------------------------
class PARAMS : public MODEL_TURB{
    public:
        #ifdef CYTHON
        PARAMS(){};// good to have for cython handling
        #endif //CYTHON
        PARAMS(string);
        //void calc_Bfield(VecDoub_I &);
        //PARAMS & operator=(const PARAMS &rhs);
};



//---------------------------------------------------
struct rhs{  
    //functor for ode; copied from Numerical Recipes
    //      Doub eps;
    //      rhs(Doub epss) : eps(epss){}
    Doub bx, by, bz;
    void operator() (PARAMS par, const Doub x, VecDoub_I &y, VecDoub_O &dydx);
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
        //void build(string, Int, Doub, Int, char*); 
        void build(const string str_tscalee, Int nsavee, int i, int j, char *dir_out);
        Output(void);
        //Output(string, const Int, char*); // mal implementado
        void init(const Int neqn, const Doub xlo, const Doub xhi, Int nHist, Int nThColl_);

        void resize(void);
        void save_dense(Stepper &, const Doub, const Doub);
        void save(const Doub, VecDoub_I &);
        void out(const Int,const Doub,VecDoub_I &,Stepper &,const Doub);
        void save2file(void);
        void claim_own(void);
        bool file_exist(void);
        #ifdef MONIT_SCATTERING
        void resizeTau(int it); // resize 'Tau[it]'
        //esto lo agrego para guardar cosas de la historia de 
        //las trayectorias:
        int nfilTau, ncolTau;       // tamanio para 'Tau'
        static const int nreg=5; // nro de regimenes de las plas
        int nreb[nreg]; // nro de rebotes/scatterings en pitch, para c/regimen
        int size_reg;
        MatDoub Tau[nreg];   // params de scattering en funcion de t
        VecDoub mu;
        Doub dtau;//collision time (instantaneous)
        void check_scattering(Doub const x, Doub *y, Doub const hdid, Doub const mu_old);
        void init_regimes(Int nHist, Int nThColl_);
        #endif // MONIT_SCATTERING
        void set_Bmodel(PARAMS*);   // para apuntar al modelo q uso en main()
        void tic(void), toc(void);  // cronometro para c/pla
        Doub trun;          // tiempo de simulacion de c/pla
        Int nsteps;         // nro total de pasos de c/pla 
        // nombres de archivos de salida 
        char fname_out[200];
        char fname_trj[200];
        char fname_misc[200];
        char fname_owned[200];

        #ifdef MONIT_STEP
        MatDoub HistStep;
        MatDoub HistSeq;
        static const Doub MaxStep=1.0;
        static const Int NStep=500;
        Doub dstep, dstep_part;
        //void monit_step(const Doub hdid);
        void monit_step(const Stepper s);
        void build_HistSeq(const Stepper s);
        #endif //MONIT_STEP

    private:
        PARAMS *pm;
        Doub vmod, bmod;
        void save_pitch(void);
        Doub bx, by, bz, vx, vy, vz;
        VecDoub XSaveGen;   // tiempos de la salida
        Int n_tscales, cc, ndec, inid, nd, base, maxd;  
        string str_tscale; //tipo de escala temporal para la salida
        Doub decade, dt;
        void set_savetimes(Doub);
        ofstream ofile;
        ofstream ofile_misc;
        ofstream ofile_own; 
    
    PRIVATE_OR_PUBLIC: // depends on CYTHON macro
        #ifdef MONIT_SCATTERING
        //----- histo del 'Tau'
        MatDoub HistTau;
        void build_HistTau(int it);
        Doub dTau, maxTau;
        Doub avrTau[nreg];
        Int nHistTau, nTau, dimHistTau;
        //int nhist_tau = nreg;
        #endif //MONIT_SCATTERING

        //--- histo del theta_coll
        Int nThColl; // has to be even!
        MatDoub HistThColl;
        void build_ThetaColl(int nr);
};


/*----- FUNCIONES NORMALES -----*/

double calc_gamma(double v);
void read_params(string fname, Doub &RIGIDITY, Doub &FRAC_GYROPERIOD, 
        Doub &NMAX_GYROPERIODS, Int &NPOINTS, Doub &ATOL, Doub &RTOL, 
        Int &n_Brealiz, string& str_timescale, Doub& tmaxHistTau, Int& nHist, Int& nThColl);
void init_orientation(int i, Doub **array_ori, VecDoub &y);
void LiberaMat(Doub **Mat, int i);
Doub **AllocMat(int nFilas, int nColumnas);
Doub **read_orientations(string fname, int &n);

#endif //FUNCS_H
//EOF
