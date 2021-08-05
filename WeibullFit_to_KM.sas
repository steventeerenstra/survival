data km;
* weibull fit to two points on survival curve;
* survival Weibull S(t)= b*exp(-(b*t)**k ) ;
t1=6; s1=0.8;
t2=18; s2=0.5;
k =( log(-log(s1)) - log(-log(s2))  )/( log(t1) - log(t2) );
b= exp( log(-log(s2))/k )/t2;
run;

** plot a curve with the above Weibull fit;
** and a curve that has a hazard that is a factor hr different;
data a;set km;
hr=1.23;
do t=1 to 30 by 1;
b2=hr*b;
S=exp(-(b*t)**k);
S2=exp( -(b2*t)**k);
output;
end;
run;

proc sgplot data=a;
series x=t y=S;
series x=t y=S2;
run;

data a; set a;
median1=(1/b)*( log(2) )**(1/k);
median2=(1/b2)*( log(2) )**(1/k);
diff_median=median1-median2;
run;

proc print data=a(obs=1); var median1 median2 diff_median;run;
*** example Baba et al. 2017;
title "Baba et al. 2017";
data km;
* weibull fit to two points on survival curve;
* survival Weibull S(t)= b*exp(-(b*t)**k ) ;
t1=6; s1=0.8;
t2=12; s2=0.5;
k =( log(-log(s1)) - log(-log(s2))  )/( log(t1) - log(t2) );
b= exp( log(-log(s2))/k )/t2;
run;

** plot a curve with the above Weibull fit;
** and a curve that has a hazard that is a factor hr different;
data a;set km;
hr=1.23;
do t=1 to 30 by 1;
b2=hr*b;
S=exp(-(b*t)**k);
S2=exp( -(b2*t)**k);
output;
end;
run;

proc sgplot data=a;
series x=t y=S;
series x=t y=S2;
run; 

data a; set a;
median1=(1/b)*( log(2) )**(1/k);
median2=(1/b2)*( log(2) )**(1/k);
diff_median=median1-median2;
run;

proc print data=a(obs=1); var median1 hr median2 diff_median;run;



*** example Kwakman et al. 2019;
title "Kwakman et al. 2019";
data km;
* weibull fit to two points on survival curve;
* survival Weibull S(t)= b*exp(-(b*t)**k ) ;
t1=6; s1=0.8;
t2=18; s2=0.5;
k =( log(-log(s1)) - log(-log(s2))  )/( log(t1) - log(t2) );
b= exp( log(-log(s2))/k )/t2;
run;

** plot a curve with the above Weibull fit;
** and a curve that has a hazard that is a factor hr different;
data a;set km;
hr=1.23;
do t=1 to 30 by 1;
b2=hr*b;
S=exp(-(b*t)**k);
S2=exp( -(b2*t)**k);
output;
end;
run;

proc sgplot data=a;
series x=t y=S;
series x=t y=S2;
run; 

data a; set a;
median1=(1/b)*( log(2) )**(1/k);
median2=(1/b2)*( log(2) )**(1/k);
diff_median=median1-median2;
run;

proc print data=a(obs=1); var median1 hr median2 diff_median;run;
