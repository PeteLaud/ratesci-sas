**********************************************************
*
* Program Name   : V_SCORECI2p.SAS
* Level / Study  : global reusable macro
* Type           : macro validation
* Description    : Validation of SCORECI macro for Poisson RD
*					(Work in progress)  
* 
* Author         : Pete Laud 
* 
* Date Created   : Dec 21 2023
* Program Status : IN DEVELOPMENT
*
**********************************************************;


*Validation of SCORECI macro, using a large random sample of parameters (alpha, n1, n0, e1, e0);
*Output results to a text/csv file for comparison against other sources:
* R 'scoreci' function for skewness-corrected interval;


*** The validation code below requires the SCORECI macro code to have been run;
*** Some environments (e.g. MS Visual Studio Code) make it difficult to get the path
*** for the current program, so the user may need to run the macro manually.
*** In interactive SAS environment, the following code should submit the macro code
*** from the directory above the location of this test program.;
%macro grabpath; 
  %qsubstr(%sysget(SAS_EXECFILEPATH),
    1,
    %length(%sysget(SAS_EXECFILEPATH)) - %length(%sysget(SAS_EXECFILENAME)) - 6
  )
%mend grabpath;
%let path = %grabpath;
%let filename = "&path.scoreci.sas";
%include &filename.;
* Otherwise, point SAS to the location of the SCORECI macro code 
* (change folder location as appropriate);
* filename prog "C:\Documents\ratesci-sas";
* %inc prog(scoreci);


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

*Run the macro on the sample of data points;
%SCORECI(DS=sample, DELTA=-0.1, LEVEL=1-alpha, STRATIFY=FALSE, SKEW=FALSE, DISTRIB=poi, OUTPUT=FALSE);
data sasval(keep=i e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set result;
 i=_n_;
run;
*Export to csv for validation against R output using program v_scoreci2p.R; 
*libname mylib "X:\My Drive\_Work\GitHub\ratesci-sas";
proc export data=sasval 
 outfile = "&path.tests\sasval2p.csv"
 dbms=csv replace;
run;

*Run again, this time with skewness correction;
%SCORECI(DS=sample, DELTA=-0.1, LEVEL=1-alpha, STRATIFY=FALSE, SKEW=TRUE, DISTRIB=poi, OUTPUT=FALSE);
*Export to csv for validation against R output using program v_scoreci.R; 
*libname mylib "X:\My Drive\_Work\GitHub\ratesci-sas";
data sasval(keep=e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set result;
 i=_n_;
run;
proc export data=sasval 
 outfile = "&path.tests\sasval2pskew.csv"
 dbms=csv replace;
run;



