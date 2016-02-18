//---
//#include "funcs.cc"
//---
template<class Stepper>
struct Odeint {
	static const Int MAXSTP=(150*50000); //MAXSTP=27*50000;
	Doub EPS;
	Int nok;
	Int nbad;
	Int nvar;			// dimension del problema: y[0], y[1], ..., y[nvar-1].
	Doub x1, x2, hmin;
	bool dense;			// True if dense output requested by out.
	VecDoub y, yold, dydx;		// (*)
	VecDoub &ystart;
	Output<Stepper> &out;
	typename Stepper::Dtype &derivs; // get the type of derivs from the stepper
	Stepper s;
	Int nstp;
	Doub x,h;
	Odeint(VecDoub_IO &ystartt,const Doub xx1,const Doub xx2,
		const Doub atol,const Doub rtol,const Doub h1,
		const Doub hminn,Output<Stepper> &outt,typename Stepper::Dtype &derivss, PARAMS, int);
	void integrate();
	PARAMS par;
	//------------------- scattering stuff
	void save_history(void);
	void check_scattering(void);
	double mu_old, mu_new, Bmod, vmod, dtau;
	// ------------------ otros
	int wrank;
// (*): en el constructor, paso las direcciones de memoria de al Stepper 's' para
// les haga las modificaciones q quiera.
};

template<class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
	const Doub atol, const Doub rtol, const Doub h1, const Doub hminn,
	Output<Stepper> &outt,typename Stepper::Dtype &derivss, PARAMS parr, int wrankk) : 
	nvar(ystartt.size()),
	y(nvar),yold(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,rtol,dense, parr), 
	par(parr), wrank(wrankk) {
	EPS=numeric_limits<Doub>::epsilon();
	h=SIGN(h1,x2-x1);
	for (Int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn, x1, x2);
}

template<class Stepper>
void Odeint<Stepper>::integrate() {
	int i=0;
	derivs(par, x, y, dydx);
	if (dense)
		out.out(-1,x,y,s,h);				// aqui solo guarda x,y
	else{
		out.save(x,y);
		i++;
		cout << " i " << i << endl;}
	dtau = 0.0;
	for (nstp=0;nstp<MAXSTP;nstp++) {
		save_history();					//--- scattering stuff
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h, derivs);

		check_scattering();				//--- scattering stuff

		if (s.hdid == h) ++nok; else ++nbad;
		if (dense){
			out.out(nstp, x, y, s, s.hdid);		// guarda solo si x>xout
		}
		else
			out.save(x,y);

		if ((x-x2)*(x2-x1) >= 0.0) {			// aqui termina la simulacion
			for (Int i=0;i<nvar;i++) ystart[i]=y[i];
			if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS){
				out.save(x,y);
			}
			out.nsteps = nstp;			// me gusta saber el nro total de pasos
			return;
		}
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
	throw("Too many steps in routine Odeint");
}

template<class Stepper>
void Odeint<Stepper>::check_scattering(){
	par.calc_Bfield(y);
	Bmod = pow(par.B[0]*par.B[0] + par.B[1]*par.B[1] + par.B[2]*par.B[2], .5);
	vmod = pow(y[1]*y[1] + y[3]*y[3] + y[5]*y[5], .5); 
	mu_new = y[1]*par.B[0] + y[3]*par.B[1] + y[5]*par.B[2];
	mu_new /= vmod*Bmod;
	//-------------------------
	dtau += s.hdid;			// controlo cuanto pasa hasta el prox rebote
	//dtauu += s.hdid;
	//-------------------------
	if(mu_old*mu_new<0.0){
		out.nreb++;
		//printf(" [rank:%d] --> nreb: %d\n", wrank, out.nreb); //getchar();
		if(out.nreb>=out.nfilTau) out.resizeTau();
		// guardo cosas de la "colisiones" con las irregularidades:
		out.Tau[out.nreb-1][0] = dtau;				// [1] tiempo-de-colision instantaneo
		out.Tau[out.nreb-1][1] = sqrt(y[0]*y[0]+y[2]*y[2]);	// [1] psic "perpend" en q ocurre dicha "colision"
		out.Tau[out.nreb-1][2] = y[4];				// [1] posic "parall"  en q ocurre dicha "colision"
		dtau = 0.0;
	}
}

template<class Stepper>
void Odeint<Stepper>::save_history(){
	/*for(int i=0; i<nvar; i++)
		yold[i] = y[i];		// guardo valores antes de integrar la ODE*/
	par.calc_Bfield(y);
	Bmod = pow(par.B[0]*par.B[0] + par.B[1]*par.B[1] + par.B[2]*par.B[2], .5);	// [G]
	vmod = pow(y[1]*y[1] + y[3]*y[3] + y[5]*y[5], .5);
	mu_old = y[1]*par.B[0] + y[3]*par.B[1] + y[5]*par.B[2];
	mu_old /= vmod*Bmod;
	//-------------------------
}

