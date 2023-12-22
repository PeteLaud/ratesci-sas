**********************************************************
*
* Program Name   : V_SCORECI.SAS
* Type           : macro validation
* Description    : Validation of SCORECI macro
*                  Including replication of published examples  
*                  and demonstration of some failings of PROC FREQ
* 
* Author         : Pete Laud 
* 
* Date Created   : Mon Mar 15 2021
* Program Status : CREATED
*
**********************************************************;

*** EXAMPLE 1
*** Demonstrating the shortcomings of standard PROC FREQ output
*** (as of SAS v9.4M7 / STAT v15.2)
*** For common risk difference, PROC FREQ gives a confidence interval
*** that is incoherent with the hypothesis test.
*** The properly stratified Miettinen-Nurminen method given by %SCORECI
*** is derived by inversion of the test statistic, so is guaranteed
*** to produce coherent results;
DATA DS1;
INPUT STRATUM   E1 N1 E0 N0;
CARDS;
 1   25 33 4  17
 2   34 50 17 25
 3   28 50 15 25
 4   52 67 22 33
;
run;
quit;

*Rearrange data for input to PROC FREQ;
data ds1t;
  set ds1;
  trt = 1; 
  outcome = 1; count = E1;  output;
  outcome = 0; count = N1 - E1; output;
  trt = 2; 
  outcome = 1; count = E0;  output;
  outcome = 0; count = N0 - E0; output;
run; 

ods output cmh=cmh 
           commonpdiff=commonpdiff 
           /*commonpdifftests=difftest*/
           crosstabfreqs=counts;
proc freq data=DS1t ;
  weight count;
  tables stratum*trt*outcome / cmh riskdiff(cl=mn common column=2) plots=riskdiffplot(CL=MN column=2 common=NO);
* Other combinations of options also do not produce the properly stratified method;
*  tables stratum*trt*outcome / cmh commonriskdiff(cl=score test=score column=2); 
*  tables stratum*trt*outcome / cmh commonriskdiff(cl=MH test=MH column=2);
run; 

title "PROC FREQ CMH test p-value";
proc print data=cmh noobs;
  var althypothesis value prob;
run;

*** Note from PROC FREQ's Method = 'Mantel-Haenszel' output, the CI (using Sato variance estimator) 
*** includes 0, contradicting the CMH test p-value (p=0.0485);
*** But with Method = 'Summary Score', the interval is shifted too far to the right
*** in comparison to the stratified MN interval (from %SCORECI macro, below);
title "PROC FREQ confidence limits for common risk difference";
proc print data=commonpdiff noobs;
  var method value lowercl uppercl;
run;

*** The next steps require the SCORECI macro code to have been run;
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

* Note the test for RD using the MN test statistic with MH weights (PVAL_2SIDED)
* is identical to the CMH test;
title "Stratified Miettinen-Nurminen confidence limits";
%SCORECI(DS=DS1, DELTA=0, LEVEL=0.95, STRATIFY=TRUE, WEIGHT=1, SKEW=FALSE);

* Applying the SKEW=TRUE option gives a skewness-corrected version of CMH, 
* reflecting a result that is slightly closer to zero than without the correction; 
title "Stratified skewness-corrected score confidence limits";
%SCORECI(DS=DS1, DELTA=0, LEVEL=0.95, STRATIFY=TRUE, WEIGHT=1, SKEW=TRUE);
title;


*** EXAMPLE 2
*** Cisapride meta-analysis dataset from Hartung & Knapp, 
*** used in Laud 2017 Appendix B (https://onlinelibrary.wiley.com/doi/10.1002/pst.1813);
*** NOTE Table BII does not show the result for the stratified CI without skewness correction;
*** (i.e. MN interval), but that can be confirmed via R ratesci package, 
*** or B.Klingenberg's code (https://sites.williams.edu/bklingen/files/2013/06/stratMHRD.r)
*** to be (0.245286, 0.369981) to 6 dps.
*** The SCAS interval to the same precision is (0.245523, 0.370330);
DATA DS2;
INPUT STRATUM   N0  E0   N1   E1;
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

*** Match published result for Method=SCAS with weights=MH;
*** (0.2456, 0.370);
*** Re-calculated in R with greater precision as
*** (0.245523, 0.370330);
%SCORECI(DS=DS2, LEVEL=0.95, STRATIFY=TRUE, SKEW=TRUE, WEIGHT=1);
%SCORECI(DS=DS2, LEVEL=0.95, STRATIFY=TRUE, SKEW=FALSE, WEIGHT=1);
%SCORECI(DS=DS2, LEVEL=0.95, STRATIFY=TRUE, SKEW=TRUE, WEIGHT=1, DISTRIB=poi);
%SCORECI(DS=DS2, LEVEL=0.95, STRATIFY=TRUE, SKEW=FALSE, WEIGHT=1, DISTRIB=poi);
proc print data=result;
 var l_bound u_bound;
run;

*** Match result from R for Method=MN with weights=MH;
*** ratesci::scoreci()
*** (0.245286, 0.369981);
%SCORECI(DS=DS2, LEVEL=0.95, STRATIFY=TRUE, SKEW=FALSE, WEIGHT=1);

*** R code for reference:
* x1 <- c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21)
* x2 <- c(9, 1, 18, 31, 6, 17, 7, 23, 3, 6, 12, 22, 19)
* n1 <- c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38)
* n2 <- c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38)
* round(scoreci(
*  x1 = x1,
*  x2 = x2,
*  n1 = n1,
*  n2 = n2,
*  stratified = TRUE, weighting = "MH", skew = TRUE, contrast = 'RD'
* )$estimates[, c(1, 3)], 6)

*Rearrange data for input to PROC FREQ, to show it gives a different CI;
data ds2t;
  set ds2;
  trt = 1; 
  outcome = 1; count = E1;  output;
  outcome = 0; count = N1 - E1; output;
  trt = 2; 
  outcome = 1; count = E0;  output;
  outcome = 0; count = N0 - E0; output;
run; 

ods output cmh=cmh 
			commonpdiff=commonpdiff 
			crosstabfreqs=counts;
proc freq data=DS2t ;
  weight count;
  tables stratum*trt*outcome / cmh riskdiff(cl=mn common column=2) plots=riskdiffplot(CL=MN column=2 common=NO);
run; 
*** Here the LCL with Method='Mantel-Haenszel' is coincidentally the same as MN to 4 dps, 
*** but the UCL is different.;

DATA DS2a;
INPUT STRATUM   N0  E0   N1   E1;
CARDS;
 1   16 1  16 12
 2  30 22 29 19
 3 29  0 56 5
 4   16 16  16 0
;
* 5 190 49 1 1 ;
*** Unstratified CIs for Laud 2017 Table BI;
%SCORECI(DS=DS2a, LEVEL=0.95, STRATIFY=FALSE, SKEW=FALSE);
%SCORECI(DS=DS2a, LEVEL=0.95, STRATIFY=FALSE, SKEW=TRUE);
%SCORECI(DS=DS2a, LEVEL=0.95, STRATIFY=FALSE, SKEW=FALSE, DISTRIB=poi);
%SCORECI(DS=DS2a, LEVEL=0.95, STRATIFY=FALSE, SKEW=TRUE, DISTRIB=poi);


*** EXAMPLE 3a: 
*** Stratified data presented in Kaifeng Lu paper 
***  (http://markstat.net/en/images/stories/lu_asa_2008.pdf);
DATA DS3a;
INPUT      STRATUM   N0  E0   N1   E1;
CARDS;
 1   15   1   5      3   
 2   10   3   10     4  
 3   25   2   35    18 
;
*** Match published result (0.1963, 0.5381); 
%SCORECI(DS=DS3a,DELTA=0.2,LEVEL=0.95,STRATIFY=TRUE, WEIGHT=1, SKEW=FALSE);


*** EXAMPLE 3b:
*** Stratified myeloma data from Klingenberg 2014;
*** https://sites.williams.edu/bklingen/files/2013/06/myel.txt;
*** Cross-check against output of Klingenberg's R code provided at
*** https://sites.williams.edu/bklingen/files/2013/06/stratMHRD.r ;
DATA DS3b;
INPUT      STRATUM  E1 N1 E0 N0;
CARDS;
1   3    4    1    3
2   3    4    8   11
3   2    2    2    3
4   2    2    2    2
5   2    2    0    3
6   1    3    2    3
7   2    2    2    3
8   1    5    4    4
9   2    2    2    3
10  0    2    2    3
11  3    3    3    3
12  2    2    0    2
13  1    4    1    5
14  2    3    2    4
15  2    4    4    6
16  4   12    3    9
17  1    2    2    3
18  3    3    1    4
19  1    4    2    3
20  0    3    0    2
21  2    4    1    5
;
*** SKEW=FALSE to demonstrate match vs published MN interval;
*** (-0.010, 0.206);
%SCORECI(DS=DS3b,DELTA=0,LEVEL=0.95,STRATIFY=TRUE, WEIGHT=1, SKEW=FALSE);
*** SKEW=TRUE to illustrate the impact of skewness correction;
*** (Note Klingenberg appears to have had problems implementing Gart-Nam method,
***  with and without skewness correction);
%SCORECI(DS=DS3b,DELTA=0,LEVEL=0.95,STRATIFY=TRUE, WEIGHT=1, SKEW=TRUE);


*** EXAMPLE 3c:
*** Stratified adverse event data from Klingenberg 2014;
*** https://sites.williams.edu/bklingen/files/2013/06/adevents.txt;
DATA DS3c;
INPUT      STRATUM  E1 N1 E0 N0;
CARDS;
1    0 432    0  142
2    1 375    0  125
3    0 80     0  40
4    0 248    1  84
5    0 50     0  24
6    2 251    0  85
7    5 605    1  209
8    1 466    2  231
9    4 525    0  175
10   1 80     0  40
11   1 375    0  124
12   0 189    0  59
13   0 55     0  27
14   0 110    0  56
15   0 85     0  85
16   0 170    0  85
17   4 1030   1  350
;
*** SKEW=FALSE to demonstrate match vs published MN interval;
*** (-0.27%, 0.38%);
%SCORECI(DS=DS3c,DELTA=0,LEVEL=0.95,STRATIFY=TRUE, WEIGHT=1, SKEW=FALSE);
*** SKEW=TRUE to illustrate the impact of skewness correction;
%SCORECI(DS=DS3c,DELTA=0,LEVEL=0.95,STRATIFY=TRUE, WEIGHT=1, SKEW=TRUE);


*** Some single stratum test cases, including examples from Newcombe 1998;
*** (https://onlinelibrary.wiley.com/doi/10.1002/%28SICI%291097-0258%2819980430%2917%3A8%3C873%3A%3AAID-SIM779%3E3.0.CO%3B2-I) ;
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
;
*** SKEW=FALSE for checking against published unstratified results from Newcombe;
%SCORECI(DS=DS4, LEVEL=0.95, delta=0, STRATIFY=FALSE, WEIGHT=1, SKEW=FALSE);

*** SKEW=TRUE to illustrate the impact of the skewness correction;
%SCORECI(DS=DS4, LEVEL=0.95, delta=0, STRATIFY=FALSE, WEIGHT=1, SKEW=TRUE);

*Rearrange data for input to PROC FREQ;
data ds4t;
  set ds4;
  trt = 1; 
  outcome = 1; count = E1;  output;
  outcome = 0; count = N1 - E1; output;
  trt = 2; 
  outcome = 1; count = E0;  output;
  outcome = 0; count = N0 - E0; output;
run; 

*Note CL=MN includes the bias correction factor N/(N-1) in the variance estimate,
* as does the CMH test,
* but the FM score test (METHOD=SCORE) does not;
ods output cmh=cmh pdifftest=pdifftest
			pdiffcls=pdiffcls 
			crosstabfreqs=counts;
proc freq data=DS4t;
  weight count;
  by stratum;
  tables trt*outcome / cmh riskdiff(cl=mn column=2 equal method=score);
run; 

*** NOTE PROC FREQ fails to produce a result for rows 5 and 6,
*** due to the data not forming a 2x2 table (0% or 100% events);
proc print data=pdiffcls;
run;


*** Trying to make sense of the one-sided tests provided by PROC FREQ;
DATA DS5;
INPUT      STRATUM   E1  N1  E0  N0;
CARDS;
 1 28 42 28 42
;
data ds5t;
  set ds5;
  trt = 1; 
  outcome = 1; count = E1;  output;
  outcome = 0; count = N1 - E1; output;
  trt = 2; 
  outcome = 1; count = E0;  output;
  outcome = 0; count = N0 - E0; output;
run; 

*** 95% MN confidence interval overlaps both -0.2 and +0.2, also in PROC FREQ;
%SCORECI(DS=DS5, LEVEL=0.95, DELTA=-0.2, STRATIFY=FALSE, WEIGHT=1, SKEW=FALSE);
ods output pdiffcls=pdiffcls cmh=cmh ;
proc freq data=DS5t;
  weight count;
  by stratum;
  tables trt*outcome / cmh riskdiff(cl=mn column=2) alpha=0.05;
run; 

*Farrington-Manning non-inferiority test p=0.0245 sig at 0.025 significance level;
*due to omission of the bias correction factor;
ods output pdiffnoninf=pdiffnoninf;
proc freq data=DS5t;
  weight count;
  by stratum;
  tables trt*outcome / cmh riskdiff(noninf margin=0.2 column=2 method=score) alpha=0.025;
run; 

*Superiority test (Pr>Z) = 0.9755 but should be (Pr<Z) = 0.0245;
ods output pdiffsup=pdiffsup;
proc freq data=DS5t ;
  weight count;
  by stratum;
  tables trt*outcome / cmh riskdiff(sup margin=0.2 column=2 method=score) alpha=0.025;
run; 

ods output pdiffequiv=pdiffequiv pdiffequivtest=pdiffequivtest;
proc freq data=DS5t;
  weight count;
  by stratum;
  tables trt*outcome / cmh riskdiff(equiv margin=0.2 column=2 method=score) alpha=0.025;
run; 
proc print data=pdiffequiv;
run;
**** Changing the margin produces a different confidence interval;
ods output pdiffequiv=pdiffequiv pdiffequivtest=pdiffequivtest;
proc freq data=DS5t;
  weight count;
  by stratum;
  tables trt*outcome / cmh riskdiff(equiv margin=0.1 column=2 method=score) alpha=0.025;
run; 
proc print data=pdiffequiv;
run;
*** It seems that the PROC FREQ "Method=Score" confidence interval uses the variance estimate evaluated at the margin,
*** which causes the CI to vary depending on the choice of MARGIN, which makes no sense;
*** CONCLUSION: Better to simply report the MN confidence interval, and conduct hypothesis tests on the same
*** underlying statistic (including bias correction) for consistency.;


*** Small sample with empty arm in a stratum;
DATA DS6;
INPUT      STRATUM   E1 N1 E0 N0;
CARDS;
 1   2 3 0 0
 2   3 5 1 2
 3   2 5 1 2
 4   5 6 2 3
;
%SCORECI(DS=DS6, DELTA=0, LEVEL=0.95, STRATIFY=TRUE, WEIGHT=1, SKEW=FALSE);
%SCORECI(DS=DS6, DELTA=0, LEVEL=0.95, STRATIFY=TRUE, WEIGHT=2, SKEW=FALSE); 


