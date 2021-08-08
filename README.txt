


1. FITTING MODEL TO DATA

$ gcc map-fitting.c -lm -I/usr/local/include -L/usr/local/lib -lgaul -lgaul_util -o map-fitting; ./map-fitting

Input files:
list_expdatafiles.dat
perturb_exp_minus10.dat
perturb_exp_minus20.dat
perturb_exp_minus30.dat
perturb_exp_minus40.dat
perturb_exp_minus50.dat
perturb_exp_plus10.dat
perturb_exp_plus20.dat
perturb_exp_plus30.dat
perturb_exp_plus40.dat
perturb_exp_plus50.dat

Output files:
hyperparams.dat
evols.dat
fitness.dat





2. SORTING AND SELECTING FITTING RESULTS

(setting N=1 selects the fittest solution; N=2 is second to fittest and so on)
$ gcc map-get-evolution.c -lm -o map-get-evolution; ./map-get-evolution N

Input file:
evols.dat
(evols_GBL.dat are the obtained solutions described in the manuscript)
fitness.dat

Output files:
evols_sort.dat
params.dat
map-plot_results.gnp
map-fitness.gnp






3. CHECKING FITTING RESULTS
3.1. RUN MODEL WITH SELECTED SOLUTION
$ gcc map-simulator.c -lm -o map-simulator; ./map-simulator

Input files:
params.dat
list_expdatafiles.dat

Output files:
perturb_model_withbl.dat
perturb_model_nobl.dat


3.2. PLOT TIME SERIES AND TRAJECTORIES
gnuplot> load 'map-plot_results.gnp'

Input file:
map-plot_results.gnp





4. RELATIONSHIP BETWEEN PARAMETER NAMES
MANUSCRIPT and MATLAB code	<--->	C CODE
a		<--->	a
b		<--->	b
c		<--->	c
d		<--->	d
alpha	<--->	delta1
beta	<--->	dseta1
gamma	<--->	eta1
delta	<--->	alpha2





0. GAUL LIBRARY
Installing gaul-devel-0.1849-0 on an Intel Core i5-7400 with Ubuntu 16.04

0.a. Tweaks before installing GAUL
At the GAUL root directory:

0.a.i. Edit file "configure"
Replace this line:
INCLUDES="$INCLUDES -I/usr/include/slang"  # FIXME: Need to detect slang.h location properly.
with this line:
INCLUDES="$INCLUDES -I/usr/include/"  # FIXME: Need to detect slang.h location properly.


0.a.ii. Edit file "/tests/test_slang.c"
Replace this line:
dief("Error %d interpreting the S-Lang script \"%s\".",SLang_Error, script_fname);
with this line:
dief("Error %d interpreting the S-Lang script \"%s\".",SLang_get_error(), script_fname);

Replace this line:
SLang_Error = 0;
with this line:
SLang_set_error(0);

Replace this line:
dief("Error %d interpreting the SLang script from stdin.", SLang_Error);
with this line:
dief("Error %d interpreting the SLang script from stdin.", SLang_get_error());

Replace this line:
SLang_Error = 0;
with this line:
SLang_set_error(0);


0.a.iii. Edit file "/util/gaul/log_util.h"
Comment out the following lines:
/*
#if defined(__GNUC__)
#  if defined(WIN32)
  FUNCPROTO inline enum log_level_type    log_get_level(void);
#  else
  FUNCPROTO extern    inline enum log_level_type  log_get_level(void);
#  endif
#else
*/

//#endif


0.a.iv. Edit file "/etc/environment"
Add the following line:
LD_LIBRARY_PATH="/usr/local/lib"


0.b. GAUL installation

./configure
make
sudo make install


