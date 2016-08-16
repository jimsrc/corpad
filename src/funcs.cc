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
    for(int i=0; i<NStep; i++){
        HistStep[i][0] = (i+.5)*dstep;
        HistStep[i][1] = 0.0;               // counts
        HistStep[i][2] = (i+.5)*(dstep_part);
        HistStep[i][3] = 0.0;               // counts
    }
    printf(" NStep: %d\n", NStep);
    printf(" MaxStep: %g\n", MaxStep);
    #endif //MONIT_STEP
}

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
        nreb++;
        if( nreb >= Tau.nrows() ) 
            resizeTau();

        //medicion de las "colisiones" con las irregularidades
        Doub* sctt = Tau[nreb-1];
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
void Output<Stepper>::init(const Int neqn, const Doub xlo, const Doub xhi, Int nHist_Tau, Int nHist_ThColl) {
    /*
     * nHist_.. : number of bins for histograms
     */
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
    }
    #ifdef MONIT_SCATTERING
    //init_regimes(nHist, nThColl_); // initilize histograms
    //-------- cosas de scatterings:
    nfilTau = 500;// tamanio inicial
    ncolTau = 5;  // 5 columnas: 1 para el tiempo, 1 para el scattering-tau, 2 para las posic parall/perp, y 1 para el angulo entre el plano x-y y z.
    Tau = MatDoub(nfilTau, ncolTau, 0.0); // (*) tiempos de scattering, y la posic x, etc.

    //-------- histograma del 'Tau'
    HistTau    = MatDoub(nHist_Tau, 2, 0.0);//(*) hist 1D
    HistThColl = MatDoub(nHist_ThColl, 2, 0.0);//(*) hist 1D
    nsteps     = 0;
    nreb       = 0;
    // (*): inicializo en ceros
    #endif // MONIT_SCATTERING
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
void Output<Stepper>::resizeTau(void){
    /* Redimensiona el vector 'xsave' hacia el doble de su longitud.
     * Como no preserva los valores, los guarda temporalmente antes de 
     * redimensionar.Despues de redimensionar,recupera la data temporal. 
     * TODO: make this an inline function! */
    nfilTau  = Tau.nrows(); // nmbr of rows NOW
    Int nold = nfilTau;
    nfilTau *= 2;
    MatDoub tempmat(Tau);   // backup de 'Tau[nr]'
    Tau.resize(nfilTau, ncolTau);
    for(Int i=0; i<nold; i++)
        for(Int j=0; j<ncolTau; j++)
            Tau[i][j] = tempmat[i][j];
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
void Output<Stepper>::build_HistTau(void){
    /* Se supone q esto se ejecuta DESPUES de q la simulacion
     * termino (todas las colisiones YA ocurrieron). Entonces
     * el 'nreb' es el nro total de colisiones ocurridas en
     * toda la simulacion. 
     * Output:
     * - HistTau : shape(nHistTau, dimHistTau=2)  
     */
    // let's find out the min/max $\tau_{coll}$ values
    Doub min=1.0e31, max=0.0, LgTau, dbin;
    Int bin;
    for(Int i=0; i<nreb; i++){
        LgTau = log10(Tau[i][1]); // log10([1/omega])
        min = MIN(min, LgTau);
        max = MAX(max, LgTau);
    }
    Int n_hist = HistTau.nrows();
    //--- build domain of histogram
    dbin = (max-min)/n_hist; //[1] define bin-width for hist
    for(int i=0; i<n_hist; i++){
        HistTau[i][0] = min + (i+.5)*dbin; // centered bins
    }
    //--- build histogram counts
    Int nxo = HistTau[0][0]/dbin; //negative index, for correction
    for(int i=0; i<nreb; i++){
        LgTau = log10(Tau[i][1]); // log10([1/omega])
        bin   = LgTau<=HistTau[0][0] ? 0 : 
                LgTau>=HistTau[n_hist-1][0] ? n_hist-1 : 
                int(LgTau / dbin) - nxo;
        HistTau[bin][1]++;
    }
    //--- calculo la media
    avrTau = 0.0;
    for(int i=0; i<nreb; i++){
        avrTau += Tau[i][1]; // [1/omega]
    }
    avrTau /= 1.0*nreb;
}


template <class Stepper>
void Output<Stepper>::build_ThetaColl(void){
    // let's find out the min/max $\theta$ values
    Doub min=1.0e31, max=0.0, theta, dbin;
    Int bin;
    // Tau[nr][:][4] ---> collision theta
    for(Int i=0; i<nreb; i++){
        theta = Tau[i][4];
        min   = MIN(min, theta);
        max   = MAX(max, theta);
    }
    //--- build domain of histogram
    dbin = (max-min)/HistThColl.nrows(); // [1] define bin-width for hist
    Int n_hist = HistThColl.nrows();
    for(int i=0; i<n_hist; i++){
        HistThColl[i][0] = min + (i+.5)*dbin; //centered bins
    }
    //--- build histogram counts
    Int nxo = HistThColl[0][0]/dbin; //negative index, for correction
    for(int i=0; i<nreb; i++){
        theta = Tau[i][4];
        bin   = theta<=HistThColl[0][0] ? 0 : 
                theta>=HistThColl[n_hist-1][0] ? n_hist-1 : 
                int(theta / dbin) - nxo;
        HistThColl[bin][1]++;
    }
    //--- calculo la media
    avrThColl = 0.0;
    for(int i=0; i<nreb; i++){
        avrThColl += Tau[i][4]; // [deg]
    }
    avrThColl /= 1.0*nreb;
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
    ofile<< "#begin_hist" << endl;
    // build histogram 'Tau[:][0,1]'
    build_HistTau(); // build histogram of $\tau_{coll}$
    //--- nro de rebotes, y colission-time promedio
    ofile<< "# n_backsctt : "<< nreb << endl; // [1]
    ofile<< "# avr_tau : "<< avrTau << endl; // [1]
    for(int i=0; i<HistTau.nrows(); i++){
        ofile << HistTau[i][0] << " ";          // [1] bin centrado
        ofile << setw(10) << HistTau[i][1] << endl; // [1] nro de cuentas en este bin
    }
    ofile << "#end_hist"<< "\n\n\n"; // (*)
    ofile<<"#END\n\n";

    //--- histograma del theta-en-colision
    ofile<<"#BEGIN THETA_COLL\n";
    ofile<<"## Histogram of angle between backscattering orientation and positive z-axis ^z+"<< endl;
    ofile<<"# avr_thcoll : "<< avrThColl << endl;
    build_ThetaColl();
    ofile<< "#begin_hist"<<endl;
    for(Int i=0; i<HistThColl.nrows(); i++){
        ofile << HistThColl[i][0] << " "; // [deg] centered bin
        ofile << setw(10) << HistThColl[i][1] << endl; // [1] counts
    }
    ofile<< "#end_hist"<< "\n\n\n";
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
