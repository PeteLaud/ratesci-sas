**********************************************************
*
* Program Name   : V_SCORECI.SAS
* Level / Study  : global reusable macro
* Type           : macro validation
* Description    : Validation of SCORECI macro
*					(Work in progress)  
* 
* Author         : Pete Laud 
* 
* Date Created   : Mon Mar 15 2021
* Program Status : IN DEVELOPMENT
*
**********************************************************;

***Validation against stratified data presented in Kaifeng Lu paper 
*  (http://markstat.net/en/images/stories/lu_asa_2008.pdf);
DATA DS1;
INPUT      STRATUM   N0  E0   N1   E1;
CARDS;
              1   15   1   5      3   
              2   10   3   10     4  
              3   25   2   35    18 
;
%SCORECI(DS=DS1,DELTA=0.2,LEVEL=0.95,STRATIFY=TRUE, WEIGHT=1, skew=FALSE);

*** Cisapride meta-analysis dataset from Hartung & Knapp, 
* used in Laud 2017 Appendix B;
DATA DS2;
INPUT      STRATUM   N0  E0   N1   E1;
CARDS;
              1   16 9  16 15
              2   16 1  16 12
              3   34 18 34 29
			  4   56 31 56 42
			  5   22 6  22 14
			  6   55 17 54 44
			  7   15 7  17 14
			  8   58 23 58 29
			  9   15 3  14 10
			  10  27 6  26 17
			  11  45 12 44 38
			  12  30 22 29 19
			  13  38 19 38 21
;
*%SCORECI(DS=DS2, LEVEL=0.95, STRATIFY=TRUE, SKEW=TRUE, WEIGHT=1);
*%SCORECI(DS=DS2, LEVEL=0.95, STRATIFY=TRUE, SKEW=TRUE, WEIGHT=2);

*** Example of incoherent results from PROC FREQ;
DATA DS3;
INPUT      STRATUM   E1 N1 E0 N0;
CARDS;
              1   25 33 4  17
              2   34 50 17 25
              3   28 50 15 25
			  4   52 67 22 33
;

*Rearrange data for input to PROC FREQ;
data ds3t;
  set ds3;
  trt=1; 
  outcome=1; count = E1;  output;
  outcome=0; count=N1-E1; output;
  trt=2; 
  outcome=1; count = E0;  output;
  outcome=0; count= N0-E0; output;
run; 

ods output cmh=cmh 
			commonpdiff=commonpdiff 
			commonpdifftests=difftest
			crosstabfreqs=counts;
proc freq data=DS3t ;
  weight count;
*  tables stratum*trt*outcome / cmh riskdiff(cl=mn common column=2);
  tables stratum*trt*outcome / cmh commonriskdiff(cl=score test=score column=2);
*  tables stratum*trt*outcome / cmh commonriskdiff(cl=MH test=MH column=2);
run; 

title "PROC FREQ CMH test p-value";
proc print data=cmh noobs;
  var althypothesis value prob;
run;

*Note with Method=Mantel-Haenszel, the CI includes 0 contradicting the test p-value;
*But with Method=Summary Score, the interval is shifted too far to the right
* in comparison to the stratified MN interval;
title "PROC FREQ confidence limits for common risk difference";
proc print data=commonpdiff noobs;
  var method value lowercl uppercl;
run;

*Note the test for RD using the MN test statistic with MH weights (PVAL_2SIDED)
* is identical to the CMH test;
title "Stratified Miettinen-Nurminen confidence limits";
%SCORECI(DS=DS3, DELTA=0, LEVEL=0.95, STRATIFY=TRUE, WEIGHT=1, skew=FALSE);
title "Stratified skewness-corrected score confidence limits";
%SCORECI(DS=DS3, DELTA=0, LEVEL=0.95, STRATIFY=TRUE, WEIGHT=1, skew=TRUE);
title;


*Some single stratum test cases, including examples from Newcombe 1998;
DATA DS4;
INPUT      STRATUM   E1  N1  E0  N0;
CARDS;
            1  56 70 48 80 
            2   9 10  3 10 
			3   6  7  2  7
			4   5 56  0 29
			5   0 10  0 20
			6  10 10 10 10
			7  10 10  0 20
			8   0 10 10 10
			9  30 75 30 30
;
%SCORECI(DS=DS4, LEVEL=0.95, STRATIFY=FALSE, WEIGHT=1, skew=FALSE, bcf=FALSE);

*Rearrange data for input to PROC FREQ;
data ds4t;
  set ds4;
  trt=1; 
  outcome=1; count = E1;  output;
  outcome=0; count=N1-E1; output;
  trt=2; 
  outcome=1; count = E0;  output;
  outcome=0; count= N0-E0; output;
run; 

*Note CL=MN includes the bias correction factor N/(N-1) in the variance estimate,
* but the score test (METHOD=SCORE) does not.;
ods output cmh=cmh 
			pdiffcls=pdiffcls pdifftest=pdifftest
			crosstabfreqs=counts;
proc freq data=DS4t ;
  weight count;
  by stratum;
  tables trt*outcome / cmh riskdiff(cl=mn equal column=2 method=score);
run; 

