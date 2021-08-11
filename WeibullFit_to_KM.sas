title "Weibull fit based on two pairs of survival time and survival probabilities (t1,s1), (t2,s2) ";
data km;
* weibull fit to two points on survival curve;
* survival Weibull S(t)= exp(-b* t**k ) ;
t1=6; s1=0.8;
t2=18; s2=0.5;
k =( log(-log(s1)) - log(-log(s2))  )/( log(t1) - log(t2) );
b= exp( log(-log(s1)) - k*log(t1) );
run;

proc print data=km; var b k;run;
** plot a curve with the above Weibull fit;
** and a curve that has a hazard that is a factor hr different;
data fit;set km;
hr=1.25;
do t=1 to 30 by 1;
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

data median; set fit;
median1=(log(2)/b)**(1/k);
median2=(log(2)/b2 )**(1/k);
diff_median=median1-median2;
relative_diff=diff_median/median1;
run;

proc print data=median(obs=1); var median1 median2 diff_median relative_diff;run;



*** example Baba et al. 2017;
title "Baba et al. 2017";
data km;
* weibull fit to two points on survival curve;
* survival Weibull S(t)= b*exp(-(b*t)**k ) ;
t1=6; s1=0.8;
t2=12; s2=0.5;
k =( log(-log(s1)) - log(-log(s2))  )/( log(t1) - log(t2) );
b= exp( log(-log(s1)) - k*log(t1) );
run;

proc print data=km; var b k;run;
** plot a curve with the above Weibull fit;
** and a curve that has a hazard that is a factor hr different;
data fit;set km;
hr=1.25;
do t=1 to 30 by 1;
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

data median; set fit;
median1=(log(2)/b)**(1/k);
median2=(log(2)/b2 )**(1/k);
diff_median=median1-median2;
relative_diff=diff_median/median1;
run;

proc print data=median(obs=1); var median1 median2 diff_median relative_diff;run;



*** example Kwakman et al. 2019;
title "Kwakman et al. 2019";
data km;
* weibull fit to two points on survival curve;
* survival Weibull S(t)= b*exp(-(b*t)**k ) ;
t1=6; s1=0.8;
t2=18; s2=0.5;
k =( log(-log(s1)) - log(-log(s2))  )/( log(t1) - log(t2) );
b= exp( log(-log(s1)) - k*log(t1) );
run;

proc print data=km; var b k;run;
** plot a curve with the above Weibull fit;
** and a curve that has a hazard that is a factor hr different;
data fit;set km;
hr=1.25;
do t=1 to 30 by 1;
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

data median; set fit;
median1=(log(2)/b)**(1/k);
median2=(log(2)/b2 )**(1/k);
diff_median=median1-median2;
relative_diff=diff_median/median1;
run;

proc print data=median(obs=1); var median1 median2 diff_median relative_diff;run;
