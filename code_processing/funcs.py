#
from pylab import *
from numpy import *
from lmfit import minimize, Parameters, Parameter, report_errors
#
def make_fit_kxx(data, sems):
        x       = data[0]
        y       = data[1]
        # create a set of Parameters
        params = Parameters()
        params.add('b')
	params.add('m')
	params.add('xo')

        SEM_b	= sems[0]
	SEM_m   = sems[1]
	SEM_xo  = sems[2]

        params['b'].value       = SEM_b
        params['b'].vary        = True
        """params['b'].min              = 
        params['b'].max         = """

        params['m'].value       = SEM_m
        params['m'].vary        = True
        """params['m'].min              = 
        params['m'].max         = """

        params['xo'].value       = SEM_xo
        params['xo'].vary        = True
        """params['xo'].min              = 
        params['xo'].max         = """

	METHOD  = "leastsq"#"leastsq"#"lbfgsb"
	result = minimize(residuals_kxx, params, args=(x, y), method=METHOD)

        # write error report
        print " --------> METODO_FITEO: %s" % METHOD
        print " --------> funcion: %s" % 'EXPRESION GIACALONE-99'
        #report_errors(params)

        par     = zeros(3)
        par[0]  = result.values['b']
	par[1]	= result.values['m']
	par[2]	= result.values['xo']
	return par

def model(N, K, L, t):
	sum = 0.
	for n in range(1, N+1):
		expon = -t * K/(L**2) * ((2.*n-1)*pi/2.)**2
		sum   += (-1.)**(n+1) / (2.*n-1) *exp(expon)
	sum *= 4./pi
	return sum

def make_fit(data, sems):
        x       = data[0]
        y       = data[1]
        # create a set of Parameters
        params = Parameters()
        params.add('N')
	params.add('A')
	params.add('yo')
        params.add('K')
        params.add('L')

        SEM_N	= sems[0]
	SEM_A   = sems[1]
	SEM_yo  = sems[2]
        SEM_K   = sems[3]
        SEM_L 	= sems[4]

        params['N'].value       = SEM_N
        params['N'].vary        = False
        """params['A'].min              = 
        params['A'].max         = """

        params['A'].value       = SEM_A
        params['A'].vary        = True
        """params['A'].min         = -.2
        params['A'].max         = +.2"""

        params['yo'].value       = SEM_yo
        params['yo'].vary        = True
	"""params['yo'].min         = -.2
        params['yo'].max         = +.2"""

        params['K'].value      = SEM_K
        params['K'].vary       = True
        """params['mu'].min     = 
        params['mu'].max        = """

        params['L'].value     = SEM_L
        params['L'].vary      = False
        """params['sig'].min    = 
        params['sig'].max       ="""

        METHOD  = "leastsq"#"leastsq"#"lbfgsb"
        result = minimize(residuals, params, args=(x, y), method=METHOD)

        # write error report
        print " --------> METODO_FITEO: %s" % METHOD
        print " --------> funcion: %s" % 'EXPRESION GIACALONE-99'
        #report_errors(params)

        par     = zeros(5)
        par[0]  = result.values['N']
	par[1]	= result.values['A']
	par[2]	= result.values['yo']
        par[3]  = result.values['K']
        par[4]  = result.values['L']
        return par

def model_shifted(N, A, yo, K, L, x):
	return yo + A*model(N, K, L, x)

def residuals(params, x, y_data):
        N      = params['N'].value
	A    = params['A'].value
	yo    = params['yo'].value
        K      = params['K'].value
        L       = params['L'].value
	diff    = (model_shifted(N, A, yo, K, L, x)  - y_data)**2.
	print " diff---> %f" % mean(diff)

        return diff

def fun_hyperbola(b, m, xo, x):
	fun = b + m/(x-xo)
	return fun

def residuals_kxx(params, x, y_data):
        b	= params['b'].value
	m	= params['m'].value
	xo	= params['xo'].value
	diff	= (fun_hyperbola(b, m, xo, x)  - y_data)**2.
	print " diff---> %f" % mean(diff)

        return diff

def load_trajectories(NBs, NPLAS, nfil, ncol, dir_data):
	data		= []
	nexist		= 0
	nwdata		= 0
	n_noexist	= 0
	for j in range(NBs):
		for i in range(NPLAS):
			fname	= '%s/B%02d_pla%03d_traj.dat' % (dir_data, j, i)
			try:
				data_aux = loadtxt(fname, unpack=True)
				nexist	+= 1
				if len(data_aux)>0:
					nwdata += 1
					data	+= [data_aux]
					print " ---> SI EXISTE: %s" % fname

			except:
				print " ---> NO EXISTE: %s" % fname
				n_noexist += 1
				#print " no existe: %s" % fname
				#pause(300)
		
	DATA 	= zeros(nwdata*nfil*ncol).reshape(nwdata*nfil, ncol)
	j	= 0
	for i in range(nwdata):
		j			= i*nfil
		DATA[j:j+nfil][:] 	= data[i].T	# eje x de 'DATA' es tiempo

	time	= data[0][0]		# [seg] con tomar una muestra es suficiente!

	# nwdata	: nro de files q existe Y tienen data
	# n_noexist	: nro de files q solicite y NO existen
	return DATA, time, nwdata, n_noexist

def sqr_deviations(DATA, time, every):
	nfil  	= len(time)
	n	= nfil/every + (nfil%every!=0)
	x2	= zeros(n)
	y2	= zeros(n)
	z2	= zeros(n)
	tt	= zeros(n)
	j	= 0
	print " calculando <x2>(t), <y2(t)>..."
	for i in range(0, nfil, every):     # yendo de 'every' en 'every'
		cond    = DATA.T[0]==time[i]
		x       = DATA.T[1][cond]
		y       = DATA.T[2][cond]
		z       = DATA.T[3][cond]
		x2[j]	= mean(x*x)	     # [AU^2]
		y2[j]	= mean(y*y)          # [AU^2]
		z2[j]	= mean(z*z)          # [AU^2]
		tt[j]	= time[i]	     # [seg]
		j 	+= 1

	return tt, x2, y2, z2



# omega ciclotron para proton:
# Ek: [eV] energ cinetica
# Bo: [G] campo guia
def calc_omega(Bo, Ek):
	q = 0.00000000048032
	mo = 1.6726e-24
	c = 3e10
	E_reposo = 938272013.0

	rigidity = sqrt(Ek*Ek + 2.*Ek*E_reposo);
	gamma = pow(pow(rigidity/E_reposo,2) + 1. , 0.5)
	beta = pow(1. - 1/(gamma*gamma) , 0.5);
	OMEGA   = q * Bo / (gamma * mo * c)		# [s^-1]
	return OMEGA

def fit_Kperp(data, fparams):
	p = make_fit_kxx(data, fparams)
	fit_b   = p[0]
	fit_m   = p[1]
	fit_xo  = p[2]
	k_fit = fit_b
	return k_fit

#
##
