# Split a multivariate ts object into univariate ts objects

`split_multi_ts()` splits a multivariate `ts` object into a list of
univariate `ts` objects, one for each variable (column).

Each resulting univariate series preserves the original time attributes
and is labeled using
[`set_ts_name()`](https://akihirao.github.io/tempssm/reference/set_ts_name.md)
with the corresponding column name. This ensures consistent handling of
variable names within the tempssm framework.

## Usage

``` r
split_multi_ts(multi_ts)
```

## Arguments

- multi_ts:

  A multivariate `ts` object. Each column represents a distinct
  variable.

## Value

A named list of univariate `ts` objects. Each list element corresponds
to one column of `multi_ts`, and the list names are taken from
`colnames(multi_ts)`.

## Details

The function requires `multi_ts` to have valid column names. If a
univariate `ts` object is supplied, the function stops with an error.

## See also

[`trim_ts_overlap`](https://akihirao.github.io/tempssm/reference/trim_ts_overlap.md),
[`set_ts_name`](https://akihirao.github.io/tempssm/reference/set_ts_name.md),
[`tempssm`](https://akihirao.github.io/tempssm/reference/tempssm.md)

## Examples

``` r
multi_ts <- ts(
  matrix(rnorm(200), ncol = 2),
  start = c(2001, 1),
  frequency = 12
)
colnames(multi_ts) <- c("precip", "solar")

split_multi_ts(multi_ts)
#> Assigned 1 variable name to time series
#> Assigned 1 variable name to time series
#> Split time series into 2 univariate series
#> $precip
#>              Jan         Feb         Mar         Apr         May         Jun
#> 2001 -0.99593086  0.42011400 -0.06286163 -0.08013185 -0.03228341 -0.71898093
#> 2002  1.51366972  0.03398894  0.74079278 -1.27080818 -0.16356423 -0.61094549
#> 2003 -0.17256490 -1.23606292 -1.90230421 -0.09450402  0.03255579  0.46129012
#> 2004 -0.88564860  0.04923707  0.18556122 -0.60865779 -0.73110285  2.71514421
#> 2005 -0.98399849  0.21472959 -0.08008508 -1.42236955  1.10814896  1.07847764
#> 2006 -1.55156017  0.77750668  1.06844014 -0.18358770  1.55824293 -0.21324238
#> 2007  1.44814671 -1.44275737  1.46718718 -0.74329990 -0.30422384  0.33765806
#> 2008  1.09173757  0.74338485 -1.20785935  0.32781805 -0.53416511  1.28394672
#> 2009 -0.04025917  1.24663211 -1.34630152 -0.57390128                        
#>              Jul         Aug         Sep         Oct         Nov         Dec
#> 2001 -1.11656132 -0.78026990 -1.77695853 -0.42783487 -2.03102701  2.75076475
#> 2002  1.25188419  0.25031947 -1.70558168 -0.85541313 -0.14490163 -0.32444696
#> 2003  1.38140030 -0.41647627  0.68094267 -0.41437304 -0.51834551 -0.68401973
#> 2004 -1.33938704 -0.64601525 -0.93245461 -0.76908693  0.37157978  0.35543278
#> 2005 -0.44025034 -0.77816901 -1.81859001 -1.12408090  1.06052384 -1.47870016
#> 2006  0.93053526  0.41081180 -1.27984430 -0.78236663 -2.27608345 -0.12639591
#> 2007 -0.60750209 -0.29556027 -0.13453714  0.81478437 -0.27292173  2.15948580
#> 2008  0.02893134 -0.39567707 -0.69486833 -1.49208070  1.44425724 -0.34701739
#> 2009                                                                        
#> 
#> $solar
#>              Jan         Feb         Mar         Apr         May         Jun
#> 2001 -0.70595925 -1.85740454 -0.24697048 -0.33556093 -0.24937687  0.45952256
#> 2002 -0.97349403 -0.51274686 -0.92348799 -1.49442573 -0.16593036  1.58223832
#> 2003 -0.19844447  0.74362012 -0.05015977  1.31314355  0.08601352 -0.62985137
#> 2004 -0.02189017 -0.83080816  1.40247432 -0.72386911  1.33058037 -0.81133302
#> 2005  1.98637658 -0.95167079 -0.55864232  0.25581454 -0.19807617  0.42878914
#> 2006 -1.13484994  0.02261431 -0.02146122 -0.27877244 -0.02420641  0.82502030
#> 2007 -0.90978818 -1.95144002 -0.80027646 -1.86936160 -0.75057148 -0.59093839
#> 2008 -0.50859339 -1.73687829  0.04048130 -0.12417946 -0.61275343  0.16066248
#> 2009  1.21387479 -0.68593246 -1.09456819  0.36487410                        
#>              Jul         Aug         Sep         Oct         Nov         Dec
#> 2001 -0.46007548  0.05893280 -0.74259108 -2.16357400 -0.23289649  0.87490683
#> 2002  1.01460227 -1.14647834 -0.90984911  3.10918268 -1.06955277 -0.70010952
#> 2003  0.65507033  0.01939355  0.69626996  2.00886989 -0.05697828  0.24041550
#> 2004  1.79503611  0.68662633  0.08937800  0.32454774  0.07131812  0.23304838
#> 2005  0.07077011 -1.53665380  0.45398437  0.87439857 -1.29773757 -1.07896524
#> 2006 -0.75018278  0.53621047 -1.57271748 -0.98642742  1.98351232 -1.85137844
#> 2007 -0.74209927  0.69351208 -0.05946397 -1.86389510 -1.27450892 -1.78177081
#> 2008 -0.66235963 -0.33485166  0.62301159  1.02819429 -1.13457824  0.91691105
#> 2009                                                                        
#> 
```
