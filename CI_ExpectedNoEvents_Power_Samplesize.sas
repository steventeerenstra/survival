** precision in HR given a number of events; 
data a; 
alpha=0.10;*two-sided (1-alpha)*100% confidence interval;
do hr=0.9, 0.95, 1.0, 1.05, 1.10;
loghr=log(hr);
do fraction= 0.2, 0.3, 0.4, 0.5, 0.6,0.8, 1.0;
total_samplesize=312+315;
events=fraction*total_samplesize;
logse=	1/sqrt(events*0.5*(1-0.5));
ul_loghr=loghr + probit(1-alpha/2)*logse;
ul_hr= exp(ul_loghr);
ll_loghr=loghr-probit(1-alpha/2)*logse;
ll_hr=exp(ll_loghr);
output;
end; end;
run;
proc print data=a noobs;by hr;var total_samplesize fraction ll_hr ul_hr;run;

** expected number of events;

data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 12 ; *accrual time;
fu=2; * followup after finishing recruitment;
median_c= 4; * median control group;
median_i= 6; * median intervention group;
N_c= 32; * sample size in the control group;
N_i= 96; * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
run;

proc print data=a;run;

* valided for OS with 15 mo vs 20 mo, 
* 32 pat/month, and 664 subjects in total, 332 vs 332,
* gives 418 deaths instead of 413 deaths in 36 total study duration, 
* so approx ok;
******************************************************************************;

** PFS of renal trial **;
data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 664/32; *accrual time: within 17 months < 664, i.e. total sample size is expected, so no followup yet;
fu=23-664/32; * followup after finishing recruitment;
median_c= 15; * median control group;
median_i= 20; * median intervention group;
N_c= 332; * sample size in the control group;
N_i=332; * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
run;

proc print data=a;run;


*********************************;
** PFS of liver trial **;
data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 30; *total accrual time, so 576/30=19.2 per month;
fu=0.5; * followup after finishing recruitment;
median_c= 8.2; * median control group;
median_i= 11.39; * median intervention group;
N_c= 192; * sample size in the control group;
N_i=2*192; * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
run;

proc print data=a;run;
 

* intermediate analysis within the accrual time;
data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 23.5 ;
fu=0; * followup after finishing recruitment;
median_c= 8.2; * median control group;
median_i= 11.39; * median intervention group;
N_c= accrual*19.2*(1/3); * sample size in the control group;
N_i=accrual*19.2*(2/3); * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
run;

proc print data=a;run;
 

** power from events, here 850 events gives 84% power for two-sided alpha is 0.025**;
data a;
alpha=0.05/2;
d=(1-0.17)*850;
theta=0.5; * fraction randomized to experimental, = 0.5 for 1:1 rando;
se=1/sqrt( d*theta*(1-theta) );
hr=0.8;
delta=abs( log(hr) );
power=probnorm(delta/se - probit(1-alpha/2) );
run;

proc print;run;

/** formulas
d= number of events;
delta=abs( log(hr) );
se= 1/sqrt(theta*(1-theta)*d)
power=probnorm(delta/se - probit(1-alpha/2) );
d=(z_(1-alpha/2) + z_(1-beta))**2  / (delta**2 * theta * (1-theta)) 
**/


*** number of events needed ***;
data a;
alpha=0.05; beta=0.05;
theta=2/3;
hr=0.5;
delta=abs(log(hr) );
d=( probit(1-alpha/2) + probit(1-beta) )**2  / ( delta**2 * theta * (1-theta)) ;
run;

proc print;run;


*** example of power calculation based on effect;
***  in change of survival rate at a landmark;
*** or in a change of medians or median and hazard ratio of course;
*input: expected survival rate in control arm S_c and intervention arm S_i ;
*       at the landmark;
*       accrual time accr, follow-up after last patient recruited fu;
*       number of patients to be recruited in ctl arm (n_ctl) and in intv arm (n_intv); 
data config;
input landmark surv_ctl_landmark surv_intv_landmark n_ctl n_intv accrual fu;
datalines;
12	0.4	0.5	80	80	36	6
12	0.4	0.55 80	80	36	6
12	0.4	0.6	70	70	36	6
;
run;

data a;set config;
alpha=0.05; *two-sided alpha;
* calculate the hazard rates in each arm and hazard ratio;
hrate_ctl= -log(surv_ctl_landmark)/landmark;
hrate_intv= -log(surv_intv_landmark)/landmark;
/** alternative in terms of change in medians;
hrate_ctl= -log(0.5)/median_ctl;
hrate_intv= -log(0.5)/median_intv;
**/
hr=hrate_intv/hrate_ctl;
median_ctl=-log(0.5)/hrate_ctl;
median_intv=-log(0.5)/hrate_intv;
*estimate the expected number of events in ctl arm at accrual+fu;
proportiondeaths_ctl=(	exp(-hrate_ctl*fu) -  exp(-hrate_ctl*(accrual+fu))	) / (hrate_ctl*accrual);
events_ctl=n_ctl*(1- proportiondeaths_ctl);
* idem in intv arm;
proportiondeaths_intv=(	exp(-hrate_intv*fu) -  exp(-hrate_intv*(accrual+fu))	) / (hrate_intv*accrual);
events_intv=n_intv*(1- proportiondeaths_intv);
* se of logrank test; 
events=events_ctl+events_intv;
theta= n_intv/(n_intv+n_ctl);
se=1/sqrt(theta*(1-theta)*events);
* power of the logrank test for the hr;
power=probnorm( abs(log(hr))/se - probit(1-alpha/2) );
run;

proc format;
value powerfmt
low -<0.70 = 'lightgray'
0.70-<0.79 = 'lightgreen'
0.79-high = 'green'
;
run;
proc tabulate data=a; var power;by landmark;
class surv_ctl_landmark surv_intv_landmark n_ctl n_intv accrual fu;
table surv_ctl_landmark*surv_intv_landmark*n_ctl*n_intv,accrual*fu*power=' '*mean=' '*[style=[background=powerfmt.]];
run;



