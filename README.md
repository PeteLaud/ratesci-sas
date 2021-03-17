# ratesci-sas

### Confidence intervals and tests for comparison of rates

ratesci-sas contains SAS macro code to compute score confidence intervals for rate (or risk) difference ('RD') for binomial proportions, with guaranteed coherence 
between the interval and the corresponding hypothesis test. [Subsequent updates may extend the code to include rate ratio ('RR', also known as relative risk), or
odds ratio ('OR'), and analysis of Poisson 'exposure-adjusted' incidence rates.]

Stratified calculations are catered for with a range of weighting schemes, with direct equivalence to the Cochran-Mantel-Haenszel (CMH) test when 
comparing RD with WEIGHT=1 for MH weights. 

Note that SAS (v9.4) PROC FREQ will produce the Miettinen-Nurminen ('MN') score interval for unstratified datasets only (and has problems producing results 
if there are no events). The "Summary Score confidence limits" produced for a stratified analysis, e.g.
 `TABLES ... / CMH COMMONRISKDIFF(CL=SCORE TEST=SCORE);`
are not stratified MN intervals, and consequently can conflict with the result of the CMH test. 
For unstratified analysis, the test and interval obtained from `TABLES ... / RISKDIFF(CL=MN EQUAL METHOD=SCORE);` 
are incoherent, because the Farrington-Manning score test omits the bias correction factor (N/(N-1)) included in the MN interval.

In addition to addressing the above issues, the SCORECI macro incorporates skewness-corrected asymptotic score ('SCAS') methods, which ensure 
improved equal-tailed coverage (or central location), 
in other words for a nominal 95% confidence interval, the one-sided non-coverage probability is (on average) close to 2.5% on each side. 
 
Corresponding hypothesis tests against any specified null parameter value (e.g. for a non-inferiority test) are provided, as well as tests of homogeneity
for stratified analysis. 

Omission of the skewness correction is also allowed, resulting in the often-recommended 'Miettinen-Nurminen' asymptotic score methods, 
which can have inferior one-sided coverage, especially for unbalanced designs. The hypothesis test for binomial RD or RR when the skewness correction is 
omitted corresponds to the Farrington-Manning test.
