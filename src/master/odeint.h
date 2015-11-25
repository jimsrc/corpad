//---
//#include "funcs.cc"
//---
template<class Stepper>
struct Odeint {
	static const Int MAXSTP=27*50000;
	Doub EPS;
	Int nok;
	Int nbad;
	Int nvar;			// dimension del problema: y[0], y[1], ..., y[nvar-1].
	Doub x1,x2,hmin;
	bool dense;			// True if dense output requested by out.
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	typename Stepper::Dtype &derivs; // get the type of derivs from the stepper
	Stepper s;
	Int nstp;
	Doub x,h;
	Odeint(VecDoub_IO &ystartt,const Doub xx1,const Doub xx2,
		const Doub atol,const Doub rtol,const Doub h1,
		const Doub hminn,Output &outt,typename Stepper::Dtype &derivss, PARAMS);
	void integrate();
	PARAMS par;
};

template<class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
	const Doub atol, const Doub rtol, const Doub h1, const Doub hminn,
	Output &outt,typename Stepper::Dtype &derivss, PARAMS parr) : nvar(ystartt.size()),
	y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,rtol,dense, parr), par(parr) {
	EPS=numeric_limits<Doub>::epsilon();
	h=SIGN(h1,x2-x1);
	for (Int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}

template<class Stepper>
void Odeint<Stepper>::integrate() {
	//printf("kmax: %d\n", out.kmax);			///
	//printf("par.b @odeint.h.integrate(): %f\n", par.b);
	int i=0;
	derivs(par,x,y,dydx);
	if (dense)
		out.out(-1,x,y,s,h);				// aqui solo guarda x,y
	else{
		out.save(x,y);
		i++;
		cout << " i " << i << endl;}
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h, derivs);
		if (s.hdid == h) ++nok; else ++nbad;
		if (dense)
			out.out(nstp, x, y, s, s.hdid);		// guarda solo si x>xout
		else
			out.save(x,y);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (Int i=0;i<nvar;i++) ystart[i]=y[i];
			if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS)
				out.save(x,y);
			return;
		}
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
	throw("Too many steps in routine Odeint");
}
