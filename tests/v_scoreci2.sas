**********************************************************
*
* Program Name   : V_SCORECI2.SAS
* Level / Study  : global reusable macro
* Type           : macro validation
* Description    : Validation of SCORECI macro
*					(Work in progress)  
* 
* Author         : Pete Laud 
* 
* Date Created   : Wed Mar 17 2021
* Program Status : IN DEVELOPMENT
*
**********************************************************;


*Validation of SCORECI macro, using a large random sample of parameters (alpha, n1, n0, e1, e0);
*Output results to a text/csv file for comparison against other sources:
*1: R 'diffscoreci' function;
*2: SAS/STAT v9.4 PROC FREQ (NB this will not calculate any interval for 0/n1-0/n0, or n1/n1-n0/n0);
*3: R 'scoreci' function for skewness-corrected interval;

%let nsamp=10000;

*Generate a random sample of numerators & denominators;
data sample(drop=treatment response freq) longsample(drop=n1 n0 e1 e0 p1 p0);
 do i=1 to &nsamp.;
  n1=ceil(rand('uniform')*300); *select n1 between 1 and 300;
  n0=ceil(rand('uniform')*300);
  e1=ceil(rand('uniform')*(n1+1))-1; *select e1 between 0 and n1;
  e0=ceil(rand('uniform')*(n0+1))-1;
  p1=e1/n1;
  p0=e0/n0;
  if mod(i,3)=0 then alpha=0.01; *select alpha level;
  else if mod(i,3)=1 then alpha=0.05;
  else alpha=0.1;
  output sample; *this wide version is the correct format for NON_INF macro to use;
  *output also in long format for proc freq to use;
   treatment="1test"; response=1; freq=e1; output longsample; *NB trts deliberately labelled in reverse because statxact does 2-1 not 1-2;
   treatment="1test"; response=2; freq=n1-e1; output longsample;
   treatment="2comp"; response=1; freq=e0; output longsample;
   treatment="2comp"; response=2; freq=n0-e0; output longsample;
 end;
run;

**Check random sample;
*proc gplot data=sample;
* plot p1*p2;
*run;
*quit;

*point SAS to the location of the NONINF macro code;
filename prog "X:\My Drive\_Work\GitHub\ratesci-sas";
%inc prog(scoreci);

* Run PROC FREQ comparison with different confidence levels;
%let alpha=0.1;
*Run the macro on the sample of data points;
%SCORECI(DS=sample,DELTA=-0.1,LEVEL=1-&alpha.,STRATIFY=FALSE, SKEW=FALSE);
data sasval(keep=i e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set result;
 i=_n_;
run;

*Check against PROC FREQ MN intervals;
ods output pdiffcls=pdiffcls;
proc freq data=longsample ;
  weight freq;
  by i;
  tables treatment*response / noprint riskdiff(cl=mn column=1) alpha=&alpha.;
run; 

data check;
 merge sasval pdiffcls (keep = i lowercl uppercl);
 by i;
 lcld = lowercl - l_bound;
 ucld = uppercl - u_bound;
run;

proc univariate data=check;
 var lcld ucld;
run;

*Run the macro again for export to R;
%SCORECI(DS=sample,DELTA=-0.1,LEVEL=1-alpha,STRATIFY=FALSE, SKEW=FALSE);
data sasval(keep=i e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set result;
 i=_n_;
run;
*Export to csv for validation against R output using program v_scoreci.R; 
*libname mylib "X:\My Drive\_Work\GitHub\ratesci-sas";
proc export data=mylib.sasval 
outfile="X:\My Drive\_Work\GitHub\ratesci-sas\sasval1.csv"
dbms=csv replace;
run;

*Run again, this time with skewness correction;
%SCORECI(DS=sample,DELTA=-0.1,LEVEL=1-alpha,STRATIFY=FALSE, SKEW=TRUE);
*Export to csv for validation against R output using program v_scoreci.R; 
*libname mylib "X:\My Drive\_Work\GitHub\ratesci-sas";
data sasval(keep=e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set result;
 i=_n_;
run;
proc export data=sasval 
outfile="X:\My Drive\_Work\GitHub\ratesci-sas\tests\sasval1skew.csv"
dbms=csv replace;
run;



