********************************************************************************;
*
* Program Name   : SCORECI.SAS
* Author         : Pete Laud, SSU, University of Sheffield
* Date Created   : 15 Mar 2021
* Version date   : 22 Dec 2023
* Repository     : Download latest version from https://github.com/PeteLaud/ratesci-sas
*
* Type           : Macro
* Description    : Computes two-sided confidence intervals (CI) for the difference
*                  between two independent binomial proportions or Poisson exposure-adjusted 
*                  rates using the method of Miettinen & Nurminen (1985), 
*                  with optional skewness correction (Laud 2017)
*                  for either stratified or unstratified 
*                  (i.e. single stratum) datasets, plus associated 
*                  two-sided superiority test (equivalent to Egon Pearson Chi-squared 
*                  test or CMH test) and/or one-sided test for a user-specified 
*                  difference, for testing non-inferiority (stratified version of 
*                  Farrington-Manning-type test).
*                  Also gives a test for homogeneity of the difference across strata 
*                  (Gart & Nam 1990).
*                  
* Macro            DS = name of input dataset, 
* variables:       DISTRIB = distribution assumed for the input data:
*                            bin = binomial (default), or poi = Poisson
*                  LEVEL = (2-sided) confidence level required, e.g. 0.95 for 95% CI
*                          NB this corresponds to a NI test at the (1-LEVEL)/2 
*                          significance level
*                  DELTA = specified non-inferiority margin (default -0.1)
*                          NB DELTA=0 corresponds to a superiority test
*                  STRATIFY = indicator for stratified or unstratified analysis 
*                           (TRUE/FALSE)
*                  SKEW   = indicator for inclusion of skewness correction for confidence interval & tests
*                           FALSE - Miettinen & Nurminen asymptotic score method
*                           TRUE (default) - Skewness-corrected asymptotic score method (SCAS)
*                  WEIGHT = weights to be used in stratified analysis:
*                           1 = MH (Mantel-Haenszel: (N1i * N0i)/(N1i + N0i), default for RD 
*                               - gives null test equivalent to CMH test, if skew=FALSE)
*                           2 = IVS (inverse variance of score  
*                               - in development, needed for Odds Ratio in a future update)
*                           3 = INV (IVS without the bias correction, 
*                               - NB this is NOT the normal approximation of 
*                                 inverse variance, but the weights from Tang 2020.
*                               - in development, for obtaining a CMH-equivalent test for OR)
*                           4 = [Minimum risk weights 
*                               - to be added, based on inverse variance of the score]
*                           5 = equal
*                           6 = user specified via WT_USER variable in dataset
*                  MAXITER, CONVERGE = precision parameters for root-finding 
*                  OUTPUT = indicator for printed output (TRUE/FALSE)
* 
* Program Status : CREATED (developed from previous NON_INF macro, 
*                           renamed appropriately for intended primary purpose)
*
* Datasets used  : Input dataset must be structured as one row per stratum<#>, 
*                  containing variables: 
*                  	N1, N0 for the sample size in the test and comparator groups
*                          (when DISTRIB = poi this is the total exposure time)
*                  	e1, e0 for the number of events in the test and 
*                              comparator groups respectively
*                  (& optional WT_USER if user-specified weights are required)
*                  <#NOTE: if stratified by more than one factor, a "stratum" is defined
*                          as a unique combination of stratification factor levels.
*                          i.e. Summing across rows should match with overall counts.>
*                  [Future update to include optional input of individual-level data.]
*						
*
* Macros used    :
*
* Output  :    RESULT dataset containing:
*                  Point estimate and confidence limits for the difference;
*                  One-sided and two-sided p-values;
*                  Maximum likelihood estimates of per-group event rates;
*              HOMTESTS dataset containing Gart-Nam test for homogeneity of stratum differences;
*              WEIGHTING dataset containing derived or provided weights;
* 
* REFERENCES:                                                                                
*        [1]. FARRINGTON, C. P. AND MANNING, G: TEST STATISTICS AND SAMPLE 
*             SIZE FORMULA FOR COMPARATIVE BINOMIAL TRIALS WITH NULL HYPOTHESIS 
*             OF NON-ZERO RISK DIFFERENCE OR NON-UNITY RELATIVE RISK,                                    
*             STATISTICS IN MEDICINE, 9:1447-1454, 1990.                                     
*                                                                                            
*        [2]. MIETTINEN, O. AND NURMINEN, M: COMPARATIVE ANALYSIS OF TWO RATES,              
*             STATISTICS IN MEDICINE, 4:213-226, 1985.                                       
*
*        [3]. Laud, P: Equal-tailed confidence intervals for comparison of rates,
*             Pharmaceutical Statistics, 16:334-348, 2017.
*
*        [4]. Tang, Y: Score confidence intervals and sample sizes for stratified 
*             comparisons of binomial proportions. 
*             Statistics in Medicine 39:3427-3457, 2020.
*
*        [5]. Gart, J and Nam, J: Approximate interval estimation of the difference in binomial 
*             parameters: Correction for skewness and extension to multiple tables. 
*             Biometrics 46(3):637-643, 1990.
*
********************************************************************************;
* 
* Amended        : Code updated to eliminate undesirable NOTEs and warnings (when options mergenoby=warn) in log;
*                : Added DF to HOMTESTS output dataset;
*                : Remove strata with N1=0 or N0=0 (for correct DF for homogeneity test);
*                : Code style improvements; 
*                : Clarified "stratum" & added repository link in header;
* Date Amended   : 22 Feb 2023
* Amended        : Added reference for Gart-Nam homogeneity test
* Date Amended   : 21 Apr 2023
* Amended        : Documentation updates in program header only
* Date Amended   : 14 Sep 2023
* Amended        : Added DISTRIB macro variable for Poisson RD calculations 
*                : Added OUTPUT macro variable to allow output to be suppressed
* Date Amended   : 22 Dec 2023
* <Repeat As Necessary following post validation amendments>
*
**********************************************************;

/********************************************************************************;
    Copyright (C) 2023 Pete Laud

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

********************************************************************************/;

OPTIONS validvarname=v7;

%MACRO SCORECI(
  DS,
  DISTRIB = bin,
  DELTA = 0,
  LEVEL = 0.95,
  STRATIFY = TRUE,
  SKEW = TRUE,
  OUTPUT = TRUE,
  WEIGHT = 1,
  MAXITER = 100,
  CONVERGE = 0.0000000001
);

* Remove any strata where the denominator is zero in either group;
data dscheck;
  set &ds.;
  if n1 = 0 or n0 = 0 then do;
    PUT "NOTE: At least one stratum contributed no information and was removed";
    delete;
  end;
run;
data dscount;
  set dscheck; 
  row = _n_;
  call symputx("nst", row); 
run;
%if "&stratify." = "TRUE"  %then %do;
  %if &nst. = 1 %then %do;
    PUT "NOTE: Only one stratum!";
	%let stratify = FALSE;
  %end;
%end;

%if "&stratify." = "FALSE" %then %do;
*for unstratified case, set weights = equal and number of strata = 1;
  %let nst = 1;
  data weights;
    set dscount(rename = (N1=N_1 N0=N_0 e1=e_1 e0=e_0));
    delta = &delta.;
	level = &level.;
    alph = 1 - level;
    W_1 = 1;
    row = _n_;
  run;
%end;

%else %if "&stratify." = "TRUE" %then %do;
 * for stratified case, get number of strata from number of rows in input dataset;
  data ds(rename = (N1=N_1 N0=N_0 e1=e_1 e0=e_0));
    attrib wt_user format = 8.4;
    set dscount; 
    row = _n_;
    if &weight. ne 6 then wt_user = 1;
    dummy = 1;
    call symputx("nst", row); 
  run;

 *convert to wide data structure;
  data wide(drop = dummy);
    merge
    %do i = 1 %to &nst.;
      ds(where = (row=&i.) 
		 rename = (N_1=N_1&i. N_0=N_0&i. e_1=e_1&i. e_0=e_0&i. WT_USER=WT_USER&i.) 
	  )
    %end;
    ;
    by dummy;
  run;

  data weights(keep = N_1: e_1: N_0: e_0: W_: delta alph level);
 				*Calculate various different weights for stratified intervals;
    set wide;
    delta = &delta.;
	level = &level.;
    alph = 1 - level;
    *set up arrays for calculations to be repeated across strata;
    array N1{&nst.} N_1:;
    array N0{&nst.} N_0:;
    array CMH{&nst.};
    array WTU(&nst.) WT_USER:;
    array W_{&nst.};

    do i = 1 to &nst.;
      CMH[i] = (N1[i] * N0[i]) / (N1[i] + N0[i]);
      if &weight. = 1 then W_[i] = CMH[i]; * /CMH_sum;   
	  	*CMH weights. NB sometimes these are called SAMPLE SIZE weights 
	  	*(e.g. Mehrotra & Railkar);
	  *[NB weight=2 & 3 are dealt with later];
      else IF (&WEIGHT. = 5 or &nst. = 1 or "&stratify." = "FALSE") THEN W_[i] = 1;  
	  			*THIS INDICATES EQUAL WEIGHTING per strata;
      else if &weight. = 6 then W_[i] = WTU[i]; *USER-SPECIFIED WEIGHTS;
    end;

  run;

  proc datasets lib = work nolist;
    delete ds wide dscheck dscount;
  run;
  quit;

%end;

*Now run the main algorithm for calculating M&N CIs;
*NB unstratified calculation is contained within this as a special case 
*with nstrat=1;
data 	main(drop = i) 
		validate(keep = N_1: N_0: e_1: e_0: W_: R1S R0S SDIFF level LL UL 
				ZZERO ZSTAT delta PVALL PVALR HOMOGQ HOMOGP incr count pt_est);
  set weights;

  *set up arrays for various calculations to be repeated across strata;
  array N1{&nst.} N_1:;
  array N0{&nst.} N_0:;
  array C1{&nst.} e_1:;
  array C0{&nst.} e_0:;
  array W_{&nst.};
  array WR1{&nst.};
  array WR0{&nst.};
  array MV{&nst.};
  array MU3{&nst.};
  array DENS{&nst.};
  array wmu3{&nst.};
  array DIFF{&nst.};

  length skew $8.;
  skew = "&skew.";

  z = probit(1 - alph/2);	*100x(LEVEL)% CUT POINT FROM Z distribution;
  pi = constant('pi'); *3.14159265358979323846264338327950288;

  *****************************************************************************;
  *WEIGHTED TEST SUBROUTINE to obtain CMH test for superiority, or 1-sided 
    "Bias-corrected stratified Farrington-Manning test" p-value
    NB The SAS/STAT v12.1 (SAS9.3) version of this test does not incorporate 
  	   the NT/(NT-1) correction term, because Farrington & Manning consider it 
  	   to be negligible for large samples. We (like STATXACT) include the 
       correction for consistency with the calculated Miettinen & Nurminen 
       confidence intervals;
  *****************************************************************************;
  D = DELTA;
  link mle;  *Solve the MLE equations as per Farrington & Manning, 
  				to obtain the denominator for the test statistic; 
  ZSTAT = ZOBS;
  PVALL = PROBNORM(ZSTAT);
  PVALR = 1 - PROBNORM(ZSTAT);
 *************************************END OF TEST SUBROUTINE********************;

  ************** CONF_INT SUBROUTINE {including point estimate} ****************;
  *Here we iterate the calculation of the score statistic, using bisection to 
  locate the 2 values of d such that ZOBS = +/-Z
  and the value of d such that ZOBS = 0;

  *initialise parameters for iteration;
  ITER = 0;

  *deal with special cases at +/- 1 not requiring iteration: 
  *e.g. if point estimate is -1 then so is LCL; 
  D = -1;
  link mle;
  Zminus1 = ZOBS;
  output main;
  D = 1;
  link mle;
  Zplus1 = ZOBS;
  output main;

  *first pass (iter = 1) finds the point estimate, then the lower CL then upper CL;
  *Finally run the point estimate again but without the skewness correction
   (for use in the homogeneity test);
  LOOP1:ITER = ITER + 1; 
  COUNT = 0;
  D1 = -1;
  D2 = 1;
  D = 0;
  if iter = 1 then Z = 0;
  else if iter = 2 then Z = probit(1 - alph/2);
  else if iter = 3 then Z = -probit(1 - alph/2);
  else if iter = 4 then do;
   Z = 0;
   skew = "FALSE";
  end;
  incr = abs(D1 - D2);
  PERR = .;
  ERR = .; 

  LOOPT:COUNT = COUNT + 1;
  IF COUNT > 1 THEN do;
  *bisection method;
    if err > 0 then D2 = D;  *If the current zobs is below the cutoff 
								then we move in the upper boundary;
    else D1 = D;           *otherwise we move in the lower boundary;
    incr = abs(D1 - D2);
    D = (D1 + D2)/2;  *Bisect the boundaries for the next iteration;
  end;

  link mle;  *Solve the MLE equations as per Farrington & Manning, 
  			to obtain the test statistic (code further down); 
  if (iter = 1 & count = 1) then Zzero = ZOBS; *save test statistic at D=0;
  PERR = ERR;  *Keep previous result for the next iteration;
  ERR = Z - ZOBS; *calculate new result and which side of the cutoff we are at;

  *deal with special cases at +/- 1 not requiring iteration: 
  *e.g. if point estimate is -1 then so is LCL; 
  if (Zminus1 <= Z) then do;
    if (iter = 1) then do; D2 = -1; incr = 0; end;
    if (iter = 2) then do; D2 = -1; incr = 0; end;
    if (iter = 4) then do; D2 = -1; incr = 0; end;
  end;
  if (Zplus1 >= Z) then do;
    if (iter = 1) then do; D1 = 1; incr = 0; end;
    if (iter = 3) then do; D1 = 1; incr = 0; end;
    if (iter = 4) then do; D1 = 1; incr = 0; end;
  end;
  if ERR = 0 then do; D1 = D; incr = 0; end;

  output main;

  IF COUNT = &maxiter. THEN 
			PUT "WARNING: CONVERGENCE NOT REACHED AFTER &maxiter. LOOPS";
  IF (incr > &converge. & COUNT < &maxiter. ) THEN do;
    GO TO LOOPT;
  end;

  *** Bisection subroutine maps [-1, +1] to [-inf, +inf] for convenience (matching code used for binomial RD);
  *** So the resulting CI for Poisson RD also needs the same transformation;
  if "&distrib." = "bin" then D1t = D1;
  else if "&distrib." = "poi" then D1t = round(tan(pi * D1 / 2), 1E-10);
  IF (ITER = 1 & (ERR = 0 | incr < &converge. | count = &maxiter.)) THEN do;
    pt_est = D1t;
    pt_back = D1;
  end;
  IF (ITER = 2 & (ERR = 0 | incr < &converge. | count = &maxiter.)) THEN LL = D1t;
  IF (ITER = 3 & (ERR = 0 | incr < &converge. | count = &maxiter.)) THEN UL = D1t;
  IF (ITER = 4 & (ERR = 0 | incr < &converge. | count = &maxiter.)) THEN uc_est = D1t;

  IF ITER < 4 THEN GOTO LOOP1;

 ********************************HOMOGENEITY TEST SUBROUTINE ****************;
  if &nst. > 1 and "&stratify." = "TRUE" then do;
    array Qstat(&nst.);
    D = pt_back; 
    link mle;
    do i = 1 to &nst.;
      Qstat[i] = ((diff[i] - pt_est)**2)/ MV[i]; *As per Laud 2017 Appendix S4.2. (S3), evaluated at MLE;
    end;
    HOMOGQ = SUM(of Qstat{*}); 
    HOMOGP = 1 - probchi(HOMOGQ, &nst.-1);
  end;
  ********************************END OF HOMOGTEST SUBROUTINE****************;

  output validate main;

  mle:
  *solve the MLE equations as per Farrington & Manning, 
  *including stratified element as per Miettinen & Nurminen;
  do i = 1 to &nst.;
    NT = N1[i] + N0[i];
    CT = C1[i] + C0[i];
    if "&distrib." = "bin" then do;
      L3 = NT;
      L2 = (N1[i] + 2*N0[i])*D - NT - CT;
      L1 = (N0[i]*D - NT - 2*C0[i])*D + CT;
      L0 = C0[i]*D*(1 - D);
      Q = (L2**3)/((3*L3)**3) - (L1*L2)/(6*(L3**2)) + L0/(2*L3);
      P = (Q >= 0)*(SQRT((L2**2)/((3*L3)**2) - L1/(3*L3)))
          - (Q < 0)*(SQRT((L2**2)/((3*L3)**2) - L1/(3*L3)));
      if Q = 0 then TEMP = 0;
      else TEMP = Q/(P**3);
      TEMP = ((1><TEMP)<>-1);  ***TO AVOID ROUNDING ERRORS IN THE ARG OF ARCOS;
      A = (1/3)*(PI + ARCOS(TEMP));
      MR0 = min(1-1E-20, max(1E-20, round(2*P*COS(A) - L2/(3*L3), 1E-12)));
      MR1 = min(1-1E-20, max(1E-20, MR0 + D));
      MV[i] = max(0, (MR1*(1 - MR1)/N1[i] + MR0*(1 - MR0)/N0[i])*(NT/(NT - 1))); 
      if (((R1S = 0 & R0S = 0) | (R1S = 1 & R0S = 1)) & D = 0) then MU3[i] = 0; 
      else MU3[i] = round((MR1*(1 - MR1)*(1 - 2*MR1)/(N1[i])**2 
                    - MR0*(1 - MR0)*(1 - 2*MR0)/(N0[i])**2), 1E-15); *Machine precision issue if e.g. MR0=0.5;
    end;
    else if "&distrib." = "poi" then do;
    * Transform D value between -1 and +1 to (-inf, +inf) for Poisson RD;
      if D = -1 then Dp = -1E12;
      else if D = 1 then Dp = 1E12;
      else Dp = round(tan(pi * D / 2), 1E-10);
      A = NT;
      B = NT*Dp - CT;
      C = -C0[i]*Dp;
      MR0 = max(1E-20, (-B + sqrt(B**2 - 4*A*C))/(2*A));
      MR1 = max(1E-20, MR0 + Dp);
      MV[i] = max(0, (MR1/N1[i] + MR0/N0[i]));
      MU3[i] = round(MR1/N1[i]**2 - MR0/N0[i]**2, 1E-15);
    end;

    if &weight. = 2 then W_[i] = 1/MV[i]; 
		* IVS weights (special handling required for zeros);
    else if &weight. = 3 then W_[i] = (1/MV[i])*(NT - 1)/NT; 
		* Alternative INV weights, Tang 2020, 
			to achieve equivalence to CMH test for OR;

  *calculate stratified point estimate;
    R1 = C1[i]/N1[i]; 
    R0 = C0[i]/N0[i];     *OBSERVED PROPORTIONS WITHIN EACH STRATUM;
    DIFF[i] = R1 - R0;      
    WR1[i] = W_[i]*R1;    *weighted observed proportions;
    WR0[i] = W_[i]*R0;
  end; 
  SUMW = sum(of W_{*}); 
  R1S = SUM(of WR1{*})/SUMW;  *THIS IS THE R1 STAR IN THE M&N PAPER (not R1 star-tilde);
  R0S = SUM(of WR0{*})/SUMW;  *THIS IS THE R0 STAR IN THE PAPER;
  SDIFF = R1S - R0S; *OBSERVED DIFF IN WEIGHTED PROPORTIONS;

  do i = 1 to &nst;
    DENS[i] = ((W_[i]/(SUMW))**2)*MV[i];
    wmu3[i] = ((W_[i]/SUMW)**3)*MU3[i];
  end;
  DEN = SUM(of DENS{*});    *DENOMINATOR FOR THE OBS CHI-SQUARE is
  							the sum of the denominators for each stratum;

  *** calculate score statistic;
  if "&distrib." = "bin" then Dt = D;
  *** Bisection subroutine maps [-1, +1] to [-inf, +inf] for convenience (matching code used for binomial RD);
  else if "&distrib." = "poi" then Dt = Dp; 
  if (SDIFF - Dt = 0) then score1 = 0;
  else score1 = (SDIFF - Dt)/sqrt(max(1E-20, DEN)); *** Avoid division by 0 NOTE because SAS cannot handle infinity;
  if (SUM(of wmu3{*}) = 0) then scterm = 0;
  else scterm = SUM(of wmu3{*})/(6*DEN**1.5);
  if (skew = "FALSE" | scterm = 0) then ZOBS = score1;
  else do;
    AA = scterm;
    BB = 1;
    CC = -1*(score1 + scterm);
    num = (-1*BB + sqrt(max(0, BB**2 - 4 * AA * CC)));
	*** For producing a warning (to be added) when a negative discriminant renders p-value non-calculable (very rare);
	dsct = BB**2 - 4 * AA * CC; 
    ZOBS = num/(2 * AA);
  end;

  return;

run;

DATA RESULT(drop = r1s r0s sdiff pt_est ll ul zstat zzero pvall pvalr homogq homogp 
					N_1: N_0: e_1: e_0:  delta level incr count w_:
				  %if "&stratify."="FALSE" %then %do; 
				    w_1 
				  %end;
				  );

  SET validate;
  DISTRIB = "&distrib.";
  SKEWCORR = "&skew.";
  CONFLEV = LEVEL;
  %if "&stratify." = "TRUE" %then %do;
    P1_AVG = R1S;
    P0_AVG = R0S;
  %end;
  %else %do;
    e1 = e_1;
    n1 = n_1;
    e0 = e_0;
    n0 = n_0;
    P1 = R1S;
    P0 = R0S;
  %end;
  TRTDIFF = PT_EST;
  L_BOUND = LL;
  U_BOUND = UL;
  CHI2 = ZZERO**2;
  PVAL_2SIDED = 1 - PROBCHI(CHI2, 1);
  test_delta = DELTA;
  Z_DIFF = ZSTAT;
  PVAL_L = PVALL;
  PVAL_R = PVALR;
RUN;

%if "&output." = "TRUE" %then %do;
  PROC PRINT DATA = RESULT noobs;
  RUN;
%end;

%if "&stratify." = "TRUE" %then %do;
DATA WEIGHTING(keep = W_: n_strata weight);
  ***TO PRESENT THE VARS IN SPECIFIC ORDER;
  LENGTH N_STRATA 8 WEIGHT $9;                      
  SET validate;
    N_STRATA = &nst.;
    IF N_STRATA = 1 THEN WEIGHT = 'NONE';            
    ELSE IF &WEIGHT. = 1 THEN WEIGHT = 'MH';         
    ELSE IF &WEIGHT. = 2 THEN WEIGHT = 'IVS';        
    ELSE IF &WEIGHT. = 3 THEN WEIGHT = 'INV';       
    ELSE IF &WEIGHT. = 5 THEN WEIGHT = 'EQUAL';      
    ELSE IF &WEIGHT. = 6 THEN WEIGHT = 'WT_USER';    
RUN;

DATA HOMTESTS(keep = Q_STAT DF PVAL_HOMOG);
  SET validate;
  Q_STAT = HOMOGQ;
  DF = &nst. - 1;
  PVAL_HOMOG = HOMOGP;
RUN;

%if "&output." = "TRUE" %then %do;
  PROC PRINT DATA = WEIGHTING noobs;
  RUN;
  PROC PRINT DATA = HOMTESTS noobs;
  RUN;
%end;

%end;

proc datasets lib = work nolist;
 delete main weights validate;
run;
quit;

%MEND SCORECI;

