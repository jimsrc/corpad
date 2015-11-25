set log x
set log y
set xlabel "time [1/omega]"
set ylabel "k_perp(t) [cm2/s]"
set grid
set macros
usin='u 1:($2+$3)'
pl 'k_vs_t_Ek.6.0e+05eV_Nm128_slab1.00_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' @usin w lp pt 7 ti "100% slab" 
repl 'k_vs_t_Ek.6.0e+05eV_Nm128_slab0.60_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' @usin w lp pt 7 ti "60% slab"
repl 'k_vs_t_Ek.6.0e+05eV_Nm128_slab0.40_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' @usin w lp pt 7 ti "40% slab"
repl 'k_vs_t_Ek.6.0e+05eV_Nm128_slab0.20_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' @usin w lp pt 4 ti "20% slab"
repl 'k_vs_t_Ek.6.0e+05eV_Nm128_slab0.10_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' @usin w lp pt 7 ti "10% slab"
repl 'k_vs_t_Ek.6.0e+05eV_Nm128_slab0.05_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' @usin w lp pt 7 ti "5% slab"
repl 'k_vs_t_Ek.6.0e+05eV_Nm128_slab0.00_sig.1.0e+00_Lc2d.3.0e-03_LcSlab.3.0e-02.dat' @usin w lp pt 7 ti "0% slab"
#    EOF
