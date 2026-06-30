# Trim and align temperature and exogenous time series over their shared period

`trim_ts_overlap()` aligns a temperature time series and one or more
exogenous time series by trimming them to their shared (overlapping)
time period.

The function returns the trimmed series as a named list of `ts` objects,
with consistent variable labels applied for downstream modeling in
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

Univariate `ts` objects are labeled using
[`set_ts_name()`](https://akihirao.github.io/tempssm/reference/set_ts_name.md)
to ensure a consistent handling of variable names within the tempssm
framework.

## Usage

``` r
trim_ts_overlap(
  temp_data,
  exo_data,
  temp_name = "temp",
  exo_name = NULL,
  temp_ts = NULL,
  exo_ts = NULL
)
```

## Arguments

- temp_data:

  A univariate `ts` object representing the temperature time series.

- exo_data:

  A multivariate or univariate `ts` object of exogenous variables. Each
  column represents a distinct exogenous covariate.

- temp_name:

  Character scalar giving the variable name for the temperature series.
  This name is applied using
  [`set_ts_name()`](https://akihirao.github.io/tempssm/reference/set_ts_name.md).

- exo_name:

  Optional character vector giving variable names for the exogenous
  variables. Its length must equal the number of exogenous variables
  when supplied. If `NULL`, default names of the form `var1`, `var2`,
  ... are assigned with a warning.

- temp_ts, exo_ts:

  Compatibility aliases for `temp_data` and `exo_data`. These are
  accepted to avoid breaking existing calls, but new code should use
  `temp_data` and `exo_data`.

## Value

A named list with the following elements:

- temperature:

  A univariate `ts` object of the trimmed temperature series.

- exogenous:

  A `ts` object of the trimmed exogenous variables.

## Details

The shared time window is determined using
[`ts.intersect()`](https://rdrr.io/r/stats/ts.union.html), and both
temperature and exogenous series are trimmed accordingly. The two inputs
must have the same frequency, but that frequency is not fixed;
quarterly, monthly, twice-monthly, and other regular `ts` frequencies
are supported.

For univariate series, variable names are assigned via
[`set_ts_name()`](https://akihirao.github.io/tempssm/reference/set_ts_name.md)
rather than [`colnames()`](https://rdrr.io/r/base/colnames.html) in
order to maintain internal consistency required by
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

## See also

[`ts.intersect()`](https://rdrr.io/r/stats/ts.union.html),
[`set_ts_name()`](https://akihirao.github.io/tempssm/reference/set_ts_name.md),
[`tempssm`](https://akihirao.github.io/tempssm/reference/tempssm.md)

## Examples

``` r
temp_data <- ts(rnorm(100), start = c(2000, 1), frequency = 12)
exo_data <- ts(matrix(rnorm(200), ncol = 2),
  start = c(2001, 1), frequency = 12
)

trim_ts_overlap(
  temp_data,
  exo_data,
  temp_name = "mean_temp",
  exo_name  = c("precip", "solar")
)
#> Assigned 1 variable name to time series
#> Trimmed series to 88 overlapping observations
#> $temperature
#>              Jan         Feb         Mar         Apr         May         Jun
#> 2001  1.22928684  0.15745303  1.40634206  0.39680092  0.32642296 -1.45934331
#> 2002 -0.54159436  0.87581342  0.74083815 -0.16481239 -0.75096428 -1.25336629
#> 2003 -1.49118975  0.81959245  1.01340453  1.05360525 -0.07142864  0.95244461
#> 2004  0.35846241  0.04818886  1.13017228 -0.52763231  0.47957533  1.41752930
#> 2005  0.59492387  1.10889090 -0.94296903  0.70034399 -0.42467483 -1.14313760
#> 2006 -0.45863225 -0.71570398  0.42954854 -1.16985274  1.32208831  1.36596464
#> 2007  2.58578811  1.21228095 -1.25119684  0.40026604  0.93470763 -0.47266985
#> 2008 -0.76936870  1.39308074  0.91866877  1.01608706                        
#>              Jul         Aug         Sep         Oct         Nov         Dec
#> 2001 -0.79929553  0.47001280 -1.15018052 -0.27159672  0.45742422 -0.01695794
#> 2002 -1.03623626 -0.02943857 -0.27494378  0.48545147  0.91698203 -0.98534392
#> 2003 -0.38403752 -0.41775999 -0.47026817 -1.69018029 -0.94549766 -0.34301196
#> 2004 -0.08009077 -0.39418700  0.52525764 -0.79593510  0.24509645 -0.77931984
#> 2005  0.23344022 -0.18923658 -1.66341126  1.91569077 -0.81689444  0.38365049
#> 2006 -0.24972371 -0.42472135  1.09843934  0.67092304  0.03095466  2.36851196
#> 2007  1.58301322  0.59775422  0.46608459 -0.74436754  0.95905526  0.15263708
#> 2008                                                                        
#> 
#> $exogenous
#>               precip        solar
#> Jan 2001  0.53463544 -0.313704285
#> Feb 2001 -0.98762856 -0.286300150
#> Mar 2001  1.18779227 -0.803440263
#> Apr 2001 -0.51754221 -0.094555455
#> May 2001 -0.25956025 -0.074046506
#> Jun 2001 -0.32806467  0.714490683
#> Jul 2001  0.07343239  0.124718366
#> Aug 2001 -0.24786302  1.296594339
#> Sep 2001 -1.37386226 -1.114492926
#> Oct 2001 -0.04044582 -0.842058813
#> Nov 2001  0.42153824 -1.504294209
#> Dec 2001  0.20159751 -0.284026455
#> Jan 2002 -1.69719192  0.042869041
#> Feb 2002  0.64228768 -0.008866413
#> Mar 2002 -0.99523961 -2.949083784
#> Apr 2002  0.96381390  0.020358254
#> May 2002 -1.65603723 -0.098446980
#> Jun 2002  1.07086109  0.589109098
#> Jul 2002 -0.10902636 -0.425751629
#> Aug 2002  1.89918639  0.634658985
#> Sep 2002 -1.13703073 -0.578524902
#> Oct 2002 -0.27971976 -0.169109048
#> Nov 2002 -0.89412905 -1.919232520
#> Dec 2002  0.13670185 -1.534266387
#> Jan 2003 -0.74916542 -1.114761221
#> Feb 2003  0.51819908  1.597811635
#> Mar 2003 -0.19233721 -0.639805078
#> Apr 2003  0.02880981  1.566690553
#> May 2003  0.35859089 -1.449155712
#> Jun 2003 -0.02899503 -0.791503959
#> Jul 2003  1.14704207 -0.504480733
#> Aug 2003  0.37358894  0.401826700
#> Sep 2003  0.32333921  0.971396469
#> Oct 2003 -0.82981932 -0.579663020
#> Nov 2003  1.39446258  1.604179789
#> Dec 2003 -0.19154358  0.225973337
#> Jan 2004  0.27227702 -0.514857122
#> Feb 2004 -1.08165901 -0.823788240
#> Mar 2004 -2.32939936  0.334415250
#> Apr 2004 -0.54962596 -0.093490765
#> May 2004 -0.07257999  0.304062078
#> Jun 2004  1.03228840 -0.476507534
#> Jul 2004  0.21513853 -0.241311717
#> Aug 2004 -0.49434012  0.824155833
#> Sep 2004  1.51213828 -1.555643970
#> Oct 2004 -0.57946326  0.093501283
#> Nov 2004  1.67478963 -0.366949421
#> Dec 2004 -1.00098026 -0.129408799
#> Jan 2005  1.22270284  0.407795520
#> Feb 2005  1.07718817 -0.583968954
#> Mar 2005 -0.61194546 -0.193449890
#> Apr 2005  0.50686697 -0.269546670
#> May 2005  0.46006837  0.073658232
#> Jun 2005  1.48439186  0.357236138
#> Jul 2005  0.88196323  0.550428418
#> Aug 2005 -0.53670224  0.038401793
#> Sep 2005  1.28555371 -1.609575292
#> Oct 2005  0.58784853 -1.049710115
#> Nov 2005 -1.30848203  2.052033974
#> Dec 2005  0.31672632  0.176436411
#> Jan 2006  1.19415830  1.128307358
#> Feb 2006  0.91305715  0.435001744
#> Mar 2006 -0.78675299  0.548833492
#> Apr 2006 -0.41092970  0.647418304
#> May 2006  0.47637368  0.878463454
#> Jun 2006  0.18591587  0.350797869
#> Jul 2006 -1.32471625  0.049879720
#> Aug 2006  1.17237106  0.835749446
#> Sep 2006  0.23149277 -0.281725199
#> Oct 2006  0.48307190 -0.791679461
#> Nov 2006 -0.53531958  0.001653727
#> Dec 2006  1.37813582 -1.187530126
#> Jan 2007 -1.30267167  0.362417465
#> Feb 2007  0.63485754 -0.549435319
#> Mar 2007  0.99965516  0.692949838
#> Apr 2007 -0.33774895 -0.060843916
#> May 2007 -0.08603361 -1.193540930
#> Jun 2007 -1.71881758 -0.119908059
#> Jul 2007 -0.92912157 -0.708107950
#> Aug 2007  0.81367605 -1.616663973
#> Sep 2007  0.52641355  0.495847824
#> Oct 2007  1.01195739  1.299750970
#> Nov 2007  0.83137981 -1.615985506
#> Dec 2007  0.41513414 -1.250367945
#> Jan 2008 -0.66125706  1.582132306
#> Feb 2008  0.54881129  1.252529118
#> Mar 2008 -0.53936780 -0.209395701
#> Apr 2008 -0.17131843 -0.639033687
#> 
```
