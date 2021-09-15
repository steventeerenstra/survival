title "Weibull fit based on two pairs of survival time and survival probabilities (t1,s1), (t2,s2) ";
data km;
* weibull fit to two points on survival curve;
* survival Weibull S(t)= exp(-b* t**k ) ;
t1=12; s1=0.553;
t2=24; s2=0.448;
k =( log(-log(s1)) - log(-log(s2))  )/( log(t1) - log(t2) );
b= exp( log(-log(s1)) - k*log(t1) );
run;

proc print data=km; var b k;run;
** plot a curve with the above Weibull fit;
** and a curve that has a hazard that is a factor hr different;
data fit;set km;
hr=0.7;
do t=0 to 60 by 1;
b2=hr*b;
S=exp(-b * t**k);
S2=exp( -b2 * t**k);
output;
end;
run;


proc sgplot data=fit;
series x=t y=S;
series x=t y=S2;
yaxis min=0 max=1;* yaxis starting at 0;
run;

data survival; set fit;
median1=(log(2)/b)**(1/k);
median2=(log(2)/b2 )**(1/k);
diff_median=median1-median2;
relative_diff=diff_median/median1;
run;
* expected medians;
proc print data=survival(obs=1); var median1 median2 diff_median relative_diff;run;


* expectedsurvival estimates;
proc print data=survival(where=( t in (6,12,18,24,36,48) )) ; var t S S2;run;

**expected number of events given a sample size;
data n_events; set survival; 
n1=300; n2=300;
e1=(1-S)*n1; e2=(1-S2)*n2;
e=e1+e2;
run;
* survival at a certain number of events;
proc print data=n_events; var n1 n2 e1 e2 e S S2; where e > 252 -10 and e < 252 +10;run;
proc print data=n_events; var n1 n2 e1 e2 e S S2; where e > 336 -10 and e < 336 +10;run;
