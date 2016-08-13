#ifndef FUNCS_CC
#define FUNCS_CC
#include "control.h"
#include "funcs.h"
#include "general.h"

using namespace std;

/*----- FUNCIONES NORMALES -----*/

/* lee parametros input en main() */
void read_params(string fname, Doub &FRAC_GYROPERIOD, 
        Doub &NMAX_GYROPERIODS, Int &NPOINTS, Doub &ATOL, Doub &RTOL, 
        Int &n_Brealiz, string& str_timescale, Doub& tmaxHistTau, 
        Int& nHist, Int& nThColl){
    string dummy;
    ifstream filein(fname.c_str());
    if (!filein.good()) {
        cout << " problema al abrir " << fname << endl;
        exit(1);
    }
    filein >> FRAC_GYROPERIOD   >> dummy;  // [1] fraction of gyroper
    filein >> NMAX_GYROPERIODS  >> dummy;  // [1] nro of gyroperiods
    filein >> NPOINTS       >> dummy;  // [1] nro output pts
    filein >> ATOL          >> dummy;  // [units of *y] abs tolerance
    filein >> RTOL          >> dummy;  // [1] rel tolerance
    filein >> n_Brealiz     >> dummy;  // [1] nmbr of B-field realizations
    filein >> str_timescale >> dummy;  // [string] nombre de la escala temporal para la salida
    filein >> tmaxHistTau   >> dummy;  // [1] max collision time (in gyro-period units) for histogram of collision times 
    filein >> nHist         >> dummy;  // [1] nmbr of bins for histogram of collision-times
    filein >> nThColl       >> dummy;  // [1] nmbr of bins for histogram of the angle between the x-y plane and z axis.
}



void init_orientation(int i, Doub **array_ori, VecDoub &y){
    double th, ph, mu;  // theta, phi y pitch-cosine
    th  = array_ori[i][0];
    ph  = array_ori[i][1];
    mu  = cos(th);  // pitch

    y[1]    = sqrt(1.-mu*mu)*cos(ph);   // [1] vx
    y[3]    = sqrt(1.-mu*mu)*sin(ph);   // [1] vy
    y[5]    = mu;               // [1] vz
    // en el origen siempre
    y[0]    = 0.0;              // [1] x
    y[2]    = 0.0;              // [1] y
    y[4]    = 0.0;              // [1] z
}



void LiberaMat(Doub **Mat, int i){ 
    for(int k=0; k<=i; k++)
            free(Mat[i]);
    free(Mat);
}



Doub **AllocMat(int nFilas, int nColumnas){
    Doub **Mat;

    if((Mat = (Doub**) calloc(nFilas, sizeof(Doub *))) == NULL)
        return NULL;

    for(int i=0; i<nFilas; i++)
        if((Mat[i] = (Doub *) calloc(nColumnas, sizeof(Doub)))==NULL) {
            LiberaMat(Mat,i);
            return NULL;
        }   
    return Mat;
}



Doub **read_orientations(string fname, int &n){
    double dummy;
    ifstream filein(fname.c_str());
    if (!filein.good()) {
        cout << " problema al abrir " << fname << endl;
        exit(1);
    }
    n = 0;
    for(;;n++){
        filein >> dummy >> dummy;
        if(filein.eof()) break;
    }

    ifstream file(fname.c_str());

    double **array;
    array   = AllocMat(n, 2);
    for(int i=0; i<n; i++){
        file >> array[i][0] >> array[i][1];
    }
    return array;
}



//----------------------- class Ouput
template <class Stepper>
void Output<Stepper>::build(const string str_tscalee, Int nsavee, int i, int j, char *dir_out){
    kmax    = 500;
    nsave   = nsavee;
    count   = 0;
    xsave.resize(kmax);
    dense   = nsave > 0 ? true : false;
    //-------------------- archivos de salida 
    sprintf(fname_out, "%s/B%02d_pla%03d", dir_out, j, i);
    sprintf(fname_trj,  "%s.dat",  fname_out);
    sprintf(fname_owned, "%s.owned", fname_out);

    str_tscale = str_tscalee;   // tipo de escala temporal para la salida

    #ifdef MONIT_STEP
    HistStep    = MatDoub(NStep, 4);
    dstep       = MaxStep/(1.0*NStep);
    dstep_part  = dstep/8.;
    //MinStep     = 1e3;
    for(int i=0; i<NStep; i++){
        HistStep[i][0] = (i+.5)*dstep;
        HistStep[i][1] = 0.0;               // counts
        HistStep[i][2] = (i+.5)*(dstep_part);
        HistStep[i][3] = 0.0;               // counts
    }
    printf(" NStep: %d\n", NStep);
    printf(" MaxStep: %g\n", MaxStep);
    //printf(" HistStep: %g\n", HistStep[0][1]);
    #endif //MONIT_STEP
}


#ifdef MONIT_SCATTERING
template <class Stepper>
void Output<Stepper>::init_regimes(Int nHist, Int nThColl_){
    /* It is assumed that XSaveGen[] is already initialized */
    //-------- cosas de scatterings:
    nfilTau = 500;// tamanio inicial
    ncolTau = 5;  // 5 columnas: 1 para el tiempo, 1 para el scattering-tau, 2 para las posic parall/perp, y 1 para el angulo entre el plano x-y y z.
    //Tau = (MatDoub*) malloc(nreg*sizeof(MatDoub));
    //TODO: REMOVE REGIMES!!! (bad idea)
    for(int i=0; i<nreg; i++){
        nreb[i] = 0; // (*) nro de rebotes en c/regimen
        Tau[i]  = MatDoub(nfilTau, ncolTau, 0.0); // (*) tiempos de scattering, y la posic x, etc... en c/ regimen
    }
    if( fmod(1.*XSaveGen.size(), 1.*nreg)!=0.0 )
        throw_nr("ERROR: nmbr of regimes must be multiple of XSaveGen.size()!");
    size_reg = XSaveGen.size()/nreg; //nmbr of XSaveGen[] covering each regime.

    //-------- histograma del 'Tau'
    nHistTau    = nHist;                   // nro de bines
    dimHistTau  = 2;
    HistTau     = MatDoub(nHistTau, dimHistTau, 0.0);//(*) histog 1-D
    nsteps      = 0;

    //-------- histograma del 'ThetaColl'
    if(fmod(nThColl_, 2.0)!=0.0) {
        printf("\n ---> ERROR: 'nThColl' has to be even!!\n");
        exit(1);
    }
    nThColl     = nThColl_;
    HistThColl  = MatDoub(nThColl, 2, 0.0); //(*) histog 1-D
    // (*): inicializo en ceros
}
#endif //MONIT_SCATTERING


#ifdef MONIT_SCATTERING
template<class Stepper>
void Output<Stepper>::check_scattering(Doub const x, Doub *y, Doub const hdid, Doub const mu_old){
    Doub xyz[3] = {y[0],y[2],y[4]};
    pm->calc_B(xyz);

    bmod = NORM(pm->B[0],pm->B[1],pm->B[2]);
    vmod = NORM(y[1],y[3],y[5]);
    Doub mu_new = (y[1]*pm->B[0]+y[3]*pm->B[1]+y[5]*pm->B[2])/(vmod*bmod);
    //-------------------------
    dtau += hdid; // controlo cuanto pasa hasta el prox rebote
    //-------------------------
    if((mu_old*mu_new)<0.0){
        int nr = int( (cc-1)/size_reg );//nr=0,1,...,nreg-1
        nreb[nr]++;
        if( nreb[nr] >= Tau[nr].nrows() ) 
            resizeTau(nr);

        // guardo cosas de la "colisiones" con las irregularidades
        Doub* sctt = Tau[nr][nreb[nr]-1];
        sctt[0] = x; //[1] time @ collision
        sctt[1] = dtau; //[1] tiempo-de-colision instantaneo
        sctt[2] = NORM(y[0],y[2],0.0); // [1] posic "perpend" en q ocurre dicha "colision"
        sctt[3] = y[4]; //[1] posic "parall" en q ocurre dicha "colision"
        sctt[4] = acos(pm->B[2]/bmod)*180./M_PI; // [deg] angulo entre plano x-y y z. Siendo 0.0 para B-vector paralelo a versor positivo ^z+.
        dtau = 0.0; // reset collision-time
    }
}
#endif //MONIT_SCATTERING


template <class Stepper>
Output<Stepper>::Output() : kmax(-1),dense(false),count(0) {}

/*
// TODO: arreglar esta implementacion si la vas a usar
template <class Stepper>
Output<Stepper>::Output(string str_tscalee, const Int nsavee, char* fname){
    build(str_tscalee, nsavee, fname);
    //kmax(500),nsave(nsavee),count(0),xsave(kmax) {
    //dense = nsave > 0 ? true : false;
}
*/

template <class Stepper>
void Output<Stepper>::set_savetimes(Doub xhi){
    if(str_tscale=="linear"){
        dxout=(x2-x1)/nsave;
        mu       = VecDoub(nsave, 0.0);     // (*)
        XSaveGen = VecDoub(nsave, 0.0);
        for(int i=0; i<nsave; i++){
            XSaveGen[i] = (i+1)*dxout;
            //printf(" XMix(%d) [wc-1]: %g\n", i, XSaveGen[i]);
        }
        cc = 0;                 // indice para 'XSaveGen'
    }
    else if(str_tscale=="mixed"){
        inid = 1;
        maxd = int(M_LOG10E*log(xhi));
        ndec = maxd - inid + 1;  // number of decades
        mu       = VecDoub(nsave*ndec, 0.0); // (*) 
        XSaveGen = VecDoub(nsave*ndec, 0.0);

        nd = inid;
        for(nd=inid; nd<maxd; nd++){
            dt = (pow(10, (1.0*(nd+1))) - pow(10, (1.0*nd)))/nsave;
            for(int i=0; i<nsave; i++){
                cc = i+(nd-inid)*nsave;
                XSaveGen[cc] = pow(10, 1.0*(nd)) + (i+1)*dt;
            }
        }

        dt = (xhi - pow(10, 1.0*(maxd) ))/nsave;
        for(int i=0; i<nsave; i++){
            cc = i+(maxd-inid)*nsave;
            XSaveGen[cc] = pow(10, 1.0*(maxd) ) + (i+1)*dt;
        }
        // reseteo indice:
        cc = 0;
    }
    else
        throw_nr(" USAR 'linear' O 'mixed' SOLAMENTE! (Jimmy)");
/*  (*): pitch en los tiempos "xsave"
 */
}

template <class Stepper>
void Output<Stepper>::init(const Int neqn, const Doub xlo, const Doub xhi, Int nHist, Int nThColl_) {
    nvar=neqn;
    if (kmax == -1) return;
    ysave.resize(nvar,kmax);
    if (dense) {
        x1=xlo;
        x2=xhi;
        xout=x1;
        if(xout>0.0) throw_nr(" NO ESTA IMPLEMENTADO PARA EMPEZAR EN TIEMPO t>0 (Jimmy).");
        //---- seteo los tiempos a guardar
        set_savetimes(xhi);
        #ifdef MONIT_SCATTERING
        init_regimes(nHist, nThColl_); // initilize histograms
        #endif // MONIT_SCATTERING
    }
}

template <class Stepper>
void Output<Stepper>::resize(){ 
    /* redimensiona el vector 'xsave' hacia el doble de su longitud, y 
       redimensiona el array 'ysave' hacia las dimensiones (nvar,kmax).
       Como no preserva los valores, los guarda temporalmente antes de 
       redimensionar.Despues de redimensionar,recupera la data temporal.*/
    Int kold=kmax;
    kmax *= 2;
    VecDoub tempvec(xsave); // backup de 'xsave'
    xsave.resize(kmax);
    for (Int k=0; k<kold; k++)
        xsave[k]=tempvec[k];
    MatDoub tempmat(ysave); // backup de 'ysave'
    ysave.resize(nvar,kmax);
    for (Int i=0; i<nvar; i++)
        for (Int k=0; k<kold; k++)
            ysave[i][k]=tempmat[i][k];
}

template <class Stepper>
void Output<Stepper>::resizeTau(int nr){
    /* Redimensiona el vector 'xsave' hacia el doble de su longitud.
     * Como no preserva los valores, los guarda temporalmente antes de 
     * redimensionar.Despues de redimensionar,recupera la data temporal. 
     * TODO: make this an inline function! */
    nfilTau  = Tau[nr].nrows(); // nmbr of rows NOW
    Int nold = nfilTau;
    nfilTau *= 2;
    MatDoub tempmat(Tau[nr]);   // backup de 'Tau[nr]'
    Tau[nr].resize(nfilTau, ncolTau);
    for(Int i=0; i<nold; i++)
        for(Int j=0; j<ncolTau; j++)
            Tau[nr][i][j] = tempmat[i][j];
}


#ifdef MONIT_STEP
template <class Stepper>
void Output<Stepper>::build_HistSeq(const Stepper s){
    HistSeq = MatDoub(s.IMAXX, 2);
    for(int i=0; i<s.IMAXX; i++){
        HistSeq[i][0] = s.nseq[i];
        HistSeq[i][1] = 0.0;
    }
}


template <class Stepper>
void Output<Stepper>::monit_step(const Stepper s){
    int ns;
    if(s.hdid<=MaxStep){
        // h total
        ns = int(s.hdid/dstep);
        HistStep[ns][1]++;
    }
    if((s.hdid/s.nstep)<=(NStep*dstep_part)){
        // h partial
        ns = int((s.hdid/s.nstep)/(dstep_part));
        HistStep[ns][3]++;
    }
    //
    for(int i=0; i<HistSeq.nrows(); i++){
        if(s.nstep==HistSeq[i][0])
            HistSeq[i][1]++;
    }
}
#endif //MONIT_STEP


template <class Stepper>
void Output<Stepper>::save_pitch(){
    Doub xyz[3] = {ysave[0][count],ysave[2][count],ysave[4][count]};
    pm->calc_B(xyz);

    bx=pm->B[0];        by=pm->B[1];        bz=pm->B[2];        // [1]
    vx=ysave[1][count]; vy=ysave[3][count]; vz=ysave[5][count]; // [1]

    bmod = NORM(bx,by,bz);
    vmod = NORM(vx,vy,vz);
    mu[count] = (vx*bx + vy*by + vz*bz)/(vmod*bmod);
}

template <class Stepper>
void Output<Stepper>::save_dense(Stepper &s, const Doub xout, const Doub h){
    if (count == kmax) resize();
    for (Int i=0;i<nvar;i++){
        ysave[i][count]=s.dense_out(i,xout,h);
    }
    save_pitch();           // calcula mu
    xsave[count++]=xout;        // <=> xsave[count]=xout; count++;
}

template <class Stepper>
void Output<Stepper>::save(const Doub x, VecDoub_I &y) {
    if (kmax <= 0) return;
    if (count == kmax) resize();
    for (Int i=0;i<nvar;i++)
        ysave[i][count]=y[i];
    save_pitch();           // calcula mu
    xsave[count++] = x;//x;     // <=> xsave[count]=x; count++;
    //cc++;
}

template <class Stepper>
void Output<Stepper>::out(const Int nstp,const Doub x,VecDoub_I &y,Stepper &s,const Doub h) {
    //if(count>=200) {printf(" COUNT=%d AQUI!\n", count); getchar();}
    if (!dense)
        throw_nr("dense output not set in Output!");
    if (nstp == -1) {
        save(x, y); 
        xout = XSaveGen[cc];
        cc++;   //+= dxout;
    } else {
        while ((x-xout)*(x2-x1) > 0.0) {
            save_dense(s, xout, h);     // interpola a 'xout'
            xout = XSaveGen[cc]; //+= dxout; // avanza 'xout' en 'dxout'
            cc++;
        }
    }
}


//esto lo agrego para guardar cosas de la historia de 
//las trayectorias:
template <class Stepper>
void Output<Stepper>::set_Bmodel(PARAMS *pmm){
    pm = pmm;
}


// escribimos archivo dummy "owned" para flagear de q YO 
// estoy trabajando con esta pla
template <class Stepper>
void Output<Stepper>::claim_own(){
    ofile_own.open(fname_owned);
    ofile_own << "dummy" << endl;
    ofile_own.close();
}


// chekea si existe (aunq su tamanio sea 0 bytes) o no un archivo.
template <class Stepper>
bool Output<Stepper>::file_exist(){
    if(ifstream(fname_owned)){
        printf("\n YA EXISTE: %s\n", fname_owned);
        return true;
    }
    printf("\n AUN NO EXISTE: %s\n", fname_owned);
    return false;
}


template <class Stepper>
void Output<Stepper>::build_HistTau(int it){
    /* Se supone q esto se ejecuta DESPUES de q la simulacion
     * termino (todas las colisiones YA ocurrieron). Entonces
     * el 'nreb' es el nro total de colisiones ocurridas en
     * toda la simulacion. 
     * Input:
     * - it: "regime" index
     * Output:
     * - HistTau : shape(nHistTau, dimHistTau=2)  */
    // reset histogram
    HistTau = MatDoub(nHistTau, dimHistTau, 0.0); // histog 1-D

    // let's find out the min/max \tau values for
    // this regime-block
    Doub min=1.0e31, max=0.0, t;
    Doub LgTau;
    for(Int i=0; i<nreb[it]; i++){
        LgTau = log10(Tau[it][i][1]); // log10([1/omega])
        min = MIN(min, LgTau);
        max = MAX(max, LgTau);
    }
    // build domain of histogram
    dTau = (max-min)/nHistTau; // [1] define bin-width for hist
    for(int i=0; i<nHistTau; i++){
        HistTau[i][0] = min + (i+.5)*dTau; // centered bins
    }
    // build histogram counts
    Int nxo = HistTau[0][0]/dTau; //negative index, for correction
    for(int i=0; i<nreb[it]; i++){
        LgTau = log10(Tau[it][i][1]); // log10([1/omega])
        nTau = LgTau<HistTau[0][0] ? 0 : 
               LgTau>HistTau[nHistTau-1][0] ? nHistTau-1 : 
               int(LgTau / dTau) - nxo;
        HistTau[nTau][1]++;
    }
    //--- calculo la media
    avrTau[it] = 0.0;
    for(int i=0; i<nreb[it]; i++){
        avrTau[it] += Tau[it][i][1]; // [1/omega]
    }
    avrTau[it] /= nreb[it];
}


template <class Stepper>
void Output<Stepper>::build_ThetaColl(int nr){
    /* Warning: 'nThColl' has to be even! 
     * NOTE: 'nr' is the id of the simulation regime  */
    HistThColl  = MatDoub(nThColl, 2, 0.0); // reset histogram
    Doub dth = 180.0/nThColl; // bin width of histogram

    // domain in th=(-90,90) [deg]
    for(int i=0; i<nThColl; i++){
        HistThColl[i][0] = -90.0 + (i+.5)*dth; // [deg]
    }

    // Tau[nr][:][4] ---> collision theta
    int nth;
    for(int i=0; i<nreb[nr]; i++){
        nth =  int(Tau[nr][i][4]/dth);
        nth += nThColl/2; // correction to avoid negative indexes
        HistThColl[nth][1]++;
    }
}


template <class Stepper>
void Output<Stepper>::save2file(){
    double t, x, y, z, v, vx, vy, vz, err;
    //-------------------- guardo la trayectoria
    ofile.open(fname_trj);
    ofile<<"#BEGIN TRAJECTORY"<<endl;
    //ofile<<"# Lc_slab : "<< pm->p_turb.Lc_slab; // [R_larmor].
    ofile<<"## format of trajectory data below:"<<endl;
    ofile<<"## t[sec]  x,y,z[AU], mu[1] (pitch), err[1] (relative error of relativistic gamma)" << endl;
    ofile<<"#begin_traj"<<endl;
    for(int i=0; i<count; i++){
        t = xsave[i];   // [1]
        x  = ysave[0][i]; y  = ysave[2][i]; z  = ysave[4][i]; // [1]
        vx = ysave[1][i]; vy = ysave[3][i]; vz = ysave[5][i]; // [1]
        v  = NORM(vx,vy,vz);  // [1]
        err = v-1.0; // velocity-relative-error (v_initial=1.0)

        ofile << setiosflags(ios :: showpoint | ios :: uppercase);
        ofile << setw(5) << setprecision(8) << t << " ";
        ofile << setw(5) << setprecision(8) << x << " ";    
        ofile << setw(5) << setprecision(8) << y << " ";   
        ofile << setw(5) << setprecision(8) << z << " ";
        ofile << setw(5) << setprecision(8) << mu[i] << " ";
        ofile << setw(5) << setprecision(8) << err << endl;
    }
    ofile<<"#end_traj"<<endl;
    // finalizamos seccion de trayectoria
    ofile<< "#END\n\n\n";

    #ifdef MONIT_SCATTERING
    /**** guardamos otras cosas sobre la historia de la trayectoria ***/
    ofile << "#BEGIN TAU_COLL\n";
    ofile << "## Histograms on measured collision-times 'Tau'" << endl;
    ofile << "# trun_minutes : "<<setw(10)<<setprecision(8)<< (trun/60.) << endl;   // [sec]
    ofile << "# steps_total : "<<setw(10)<<setprecision(10)<< nsteps << endl; // total nmbr of steps
    ofile << "# n_regimes : "<< nreg << endl; //nmbr of regimes
    //TODO: REMOVE REGIMES!! (bad idea)
    for(int nr=0; nr<nreg; nr++){
        //--- histogramas de 'Tau[nr][:][1]' 
        ofile<< "#begin_hist_"<< nr << endl;
        build_HistTau(nr); // build histogram for regime 'nr'
        //--- nro de rebotes, y colission-time promedio
        ofile<< "# n_rebotes : "<< nreb[nr] << endl; // [1]
        ofile<< "# tau_avr : "<< avrTau[nr] << endl; // [1]
        for(int i=0; i<HistTau.nrows(); i++){
            ofile << HistTau[i][0] << " ";          // [1] bin centrado
            ofile << setw(10) << HistTau[i][1] << endl; // [1] nro de cuentas en este bin
        }
        ofile << "#end_hist_"<< nr << "\n\n\n"; // (*)
    }
    ofile<<"#END\n\n";

    //--- histograma del theta-en-colision
    ofile<< "#BEGIN THETA_COLL\n";
    ofile<< "## Histogram of angle between x-y plane and z axis (Theta_Coll)"<< endl;
    for(int nr=0; nr<nreg; nr++){
        build_ThetaColl(nr);
        ofile<< "#begin_hist_ThetaColl_"<< nr <<endl;
        for(Int i=0; i<nThColl; i++){
            ofile << HistThColl[i][0] << " "; // [deg] bin centrado
            ofile << setw(10) << HistThColl[i][1] << endl; // [1] nro de cuentas
        }
        ofile<< "#end_hist_ThetaColl_"<< nr << "\n\n\n";
    }
    ofile<<"#END\n";
    #else
    ofile << "## +++++ NO SCATTERING INFORMATION +++++" << endl;
    #endif //MONIT_SCATTERING

    // cerramos archivo de misc
    ofile.close();
    /* (*): to make it gnuplot-friendly  */
}


template <class Stepper>
void Output<Stepper>::tic(){
    trun = time(NULL);
}

template <class Stepper>
void Output<Stepper>::toc(){
    trun = time(NULL) - trun; // [sec] tiempo de corrida para 1 pla
}



//------------------------------------------- class PARAMS
PARAMS::PARAMS(string fname_turb):
    MODEL_TURB(fname_turb) {
}



//-------------------------------------------
void rhs::operator() (PARAMS par, const Doub x, VecDoub_I &y, VecDoub_O &dydx ){
    //double bx, by, bz; 
    Doub xyz[3] = {y[0],y[2],y[4]};
    par.calc_B(xyz);
    bx=par.B[0]; by=par.B[1]; bz=par.B[2];
    // rewrite x^2y"(x)+xy'(x)+x^2y=0 as coupled FOODEa
    dydx[0] = y[1];
    dydx[1] = y[3] * bz - y[5] * by; 
    dydx[2] = y[3];
    dydx[3] =-y[1] * bz + y[5] * bx; 
    dydx[4] = y[5];
    dydx[5] =-y[3] * bx + y[1] * by; 
}


// declare/define (?) class w/ specific template
#include "stepperbs.h"
template class Output<StepperBS<rhs> >; // rhs: system of equations I use!


#endif //FUNCS_CC
//EOF
