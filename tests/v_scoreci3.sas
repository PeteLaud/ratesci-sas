**********************************************************
*
* Program Name   : V_SCORECI3.SAS
* Level / Study  : global reusable macro
* Type           : macro validation
* Description    : Validation of SCORECI macro with STRATIFY=TRUE
*					(Work in progress)  
* 
* Author         : Pete Laud 
* 
* Date Created   : Oct 29 2023
* Program Status : IN DEVELOPMENT
*
**********************************************************;

*** Validation of SCORECI macro, using a random sample of parameters (alpha, n1, n0, e1, e0);
*** Output results to a text/csv file for comparison against R 'scoreci' function.
***  NOTE SAS/STAT PROC FREQ currently (as of SAS/STAT 15.2 2023) 
***  does *not* provide the stratified M-N method;

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

*** Select number of samples, number of strata, and stratum size;
%let nsamp = 100;
%let nstrat = 4;
%let maxn = 10;

*Generate a random sample of numerators & denominators;
data samples; 
 do sample = 1 to &nsamp.;
  do stratum = 1 to &nstrat.;
   n1 = ceil(rand('uniform')*&maxn.); *select n1 between 1 and &maxn;
   n0 = ceil(rand('uniform')*&maxn.);
   e1 = ceil(rand('uniform')*(n1+1))-1; *select e1 between 0 and n1;
   e0 = ceil(rand('uniform')*(n0+1))-1;
   p1 = e1/n1;
   p0 = e0/n0;
   if mod(sample,3)=0 then alpha=0.01; *select alpha level;
   else if mod(sample,3)=1 then alpha=0.05;
   else alpha=0.1;
   output samples; *this wide version is the correct format for SCORECI macro to use;
  end;
 end;
run;



%macro runstrat(skew=FALSE);

%do sample = 1 %to &nsamp.;

 data onesample;
  set samples;
  if sample = &sample.;
  call symput('alpha', alpha);
 run;

 %SCORECI(DS=onesample, DELTA=-0.1, LEVEL=1-&alpha., STRATIFY=TRUE, SKEW=&skew.);
 
 data oneresult;
  sample = &sample.;
  set result;
 run; 

 proc append base=results data=oneresult force;
 run;
 quit;

 proc datasets lib=work nolist;
  delete homtests weighting result oneresult onesample;
 run;
 quit;

%end;
%mend;

%runstrat(skew = FALSE);

data output;
 merge samples results;
 by sample;
run;

*Export to csv for validation against R output using program v_scoreci3.R; 
data sasval(keep=sample stratum e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set output;
run;

*Export to csv for validation against R output using program v_scoreci.R; 
proc export data=sasval 
 outfile = "&path.tests\sasval3.csv"
 dbms = csv replace;
run;

 proc datasets lib=work nolist;
  delete output results sasval;
 run;
 quit;

* (Repeat with SKEW=TRUE);
%runstrat(skew=TRUE);

data output;
 merge samples results;
 by sample;
run;

*Export to csv for validation against R output using program v_scoreci3.R; 
data sasval(keep=sample stratum e1 n1 e0 n0 conflev l_bound u_bound test_delta pval_L);
 set output;
run;

*Export to csv for validation against R output using program v_scoreci.R; 
*libname mylib "X:\My Drive\_Work\GitHub\ratesci-sas";
proc export data=sasval 
 outfile = "&path.tests\sasval3skew.csv"
 dbms = csv replace;
run;

 proc datasets lib=work nolist;
  delete output results sasval;
 run;
 quit;
