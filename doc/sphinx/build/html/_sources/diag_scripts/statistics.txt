:mod:`statistics`
=================
.. function:: dim_stddev_wgt_Wrap(field[*]:numeric, ww[*]:numeric, opt[1]:integer)

   :param numeric field: a one-dimensional numeric array.
   :param numeric ww: a one-dimensional numeric array of the same size of field.
   :param integer opt: a scalar, it has the same meaning as in the corresponding NCL function dim_avg_wgt_Wrap

   Return value
      A float or a double depending on the type of input
  
   Description
      Calculates the (unbiased) weighted standard deviation consistently with
      the NCL function dim_std_dev (i.e. it divides by N-1 instead of N). For
      the weighted case this means applying a correction factor:
        sum(w_i)/[sum(w_i)^2 - sum(w_i^2)]
      Missing values are ignored.
  
   Caveats
  
   References
      en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
  
   Modification history
      20141215-A_righ_ma: written.
  
.. function:: time_operations(field:numeric, y1[1]:integer, y2[1]:integer, oper[1]:string, opt[1]:string, l_wgt[1]:logical)

   :param numeric field: a numeric array of rank 1 to 4, first dimension must be time.
   :param integer y1: start year of the time period to be averaged (-1 for full range).
   :param integer y2: end year of the time period to be averaged (-1 for full range).
   :param string oper: type of operations: "extract": no average, just extract selected period. "average": average. "stddev": (unbiased) standard deviation.
   :param string opt: operation options (has no effect is oper = extract): "annualclim": annual climatology. "seasonalclim": seasonal climatology for the standard seasons DJF, MAM, JJA, SON. "monthlyclim": monthly climatology jan-dec. For monthly input data only! Apply mymm first, if necessary. "mymm": multi year monthly mean "yearly": time average over every year in [y1:y2]. [month strings]: climatology of selected (consecutive) months (e.g., "MAM", "SONDJ"). [1, 12]: climatology of the selected month ("1"=Jan, "2"=Feb, ..., "12"=Dec).
   :param logical l_wgt: if True, calculate weighted average, with days-per-month as weights (has no effect is opt = "extract").

   Return value
      An array of the same rank as field or of rank-1, depending on oper/opt.
  
   Description
      Performs differnt types of time average, standard deviation or extraction
      of a selected time period. Weighted average (with days-per-month as
      weights) can be optionally applied.
  
   Caveats
      The weighted standard deviation is not yet implmented for all cases
      The weighted standard deviation is calculated using the unbiased
      estimator, c
      This should take into account missing values and exclude the w_i for
      which the field contains only missing values. This feature is not
      implemented yet.
  
   References
  
   Modification history
      20140703-A_gott_kl: added option "mymm".
      20140312-A_righ_ma: extended with standard deviation.
      20140109-A_righ_ma: written.
  
.. function:: calc_season_index(season[1]:string)

   :param string season: the season in upper case.

   Return value
      The indices to the months in season, e.g. "JFM" returns (/0, 1, 2/).
  
   Description
      Given the "season", i.e., any substring from "JFMAMJJASOND", retrieves
      the corresponding indices. Crashes if given substring is not unique or
      does not exist.
  
   Caveats
  
   References
  
   Modification history
  
.. function:: extract_season(data:numeric, season[1]:string)

   :param numeric data: a numeric field with time dimension.
   :param string season:  the season in upper case.

   Return value
      The temporal subset of indata defined by the 'season' string.
  
   Description
      Given the "season", i.e., any substring from "JFMAMJJASOND", retrieves
      the corresponding months from data.
  
   Caveats
  
   References
  
   Modification history
  
.. function:: month_to_season_extended(indata[*][*][*]:float, season[1]:string)

   :param float indata: a [lat][lon][time] array.
   :param string season: compute the average for this season.

   Return value
      An array with the seasonal average for each year.
  
   Description
      For each year in the input data, averages indata over the given season.
  
   Caveats
  
   References
  
   Modification history
  
.. function:: coswgt_areaave(field:numeric)

   :param numeric field: numeric field.

   Return value
      The area average using cosine lat weights.
  
   Description
      Computes the area average using cosine lat weights and lon weights=1.
  
   Caveats
  
   References
  
   Modification history
      20131209-A_eval_ma: written.
  
.. function:: coswgt_arearmse(field1:numeric, field2:numeric)

   :param numeric field1: numeric field
   :param numeric field2: numeric field

   Return value
      Area rmse average using cosine lat weights.
  
   Description
      Computes area rmse areage using cosine lat weights and lon weights=1.
  
   Caveats
  
   References
  
   Modification history:
      20131209-A_eval_ma: written.
  
.. function:: coswgt_pattern_cor(field1:numeric, field2:numeric)

   :param numeric field1: numeric field.
   :param numeric field2: numeric field.

   Return value
      Pattern correlation cosine lat weights.
  
   Description
  
   Caveats
  
   References
  
   Modification history:
      20140115-A_eval_ma: written.
  
.. function:: interannual_variability(field: numeric, y1[1]: integer, y2[1]: integer, opt[1]: string)

   :param  numeric field: a numeric array of rank 1 to 4, first dimension must be time.
   :param  integer y1: start year of the time period to be averaged (-1 for full range).
   :param  integer y2: end year of the time period to be averaged (-1 for full range).
   :param  string opt: operation options (same as time_operations): "annualclim": annual climatology. "seasonalclim": seasonal climatology for the standard seasons DJF, MAM, JJA, SON. "monthlyclim": monthly climatology jan-dec. [month strings]: climatology of selected (consecutive) months (e.g. "MAM", "SONDJ"). [1, 12]: climatology of the selected month ("1"=Jan, "2"=Feb, ..., "12"=Dec).

   Return value
      An array of the same rank as field or of rank-1, depending on opt.
  
   Description
      Calculates the standard deviation with respect to interannual
      variability, to be used as input for statistical tests.
  
   Caveats
      The standard deviation is not weighted, being w.r.t. interannual
      variability for which all years have the same weight.
  
   Reference
  
   Modification history
      20140314-A_righ_ma: written.
  
.. function:: calculate_metric(var:numeric, ref:numeric, metric:string)

   :param numeric var: a 1-D or 2-D numerical array.
   :param numeric ref: a numerical array of the same dimensionality of var.
   :param string metric: a string with the metric to calculate: "RMSD": root-mean square difference. "BIAS": mean bias. "stddev_ratio": ratio of standard deviations of var and ref (to be used in Taylor diagram). "correlation": pattern correlation for var and ref (to be used in Taylor diagram).

   Return value
      A scalar float representing the calculated grading metric.
  
   Description
      Calculate a grading metrics given two input variables of the same
      dimensionality.
  
   Modification history
      20140313-A_righ_ma: implemented weights calculation within the function,
                       depending on dimensionality.
      20140120-A_fran_fr: written.
  
.. function:: normalize_metric(var:numeric, opt:string)

   :param numeric var: numerical array.
   :param string opt: option determining the used normalization: "mean": normalization with mean. "median": normalization with median. "stddev_mean": normalization with substracting the mean and dividing by the standard deviation. "centered_median": substracting and dividing by the median.

  
   Return value
      A numerical array of the same dimensionality as var.
  
   Description
      Normalizes an array of metrics according to opt.
  
   Caveats
      Treatment of missing values not explicitely specified (yet).
  
   Reference
  
   Modification history
      20140609-A_righ_ma: absolute value added to "mean" normalization.
      20140120-A_fran_fr: written.
  
.. function:: distrib_stats(var[*]:numeric, opt:string)

   :param numeric var: a one-dimensional input array.
   :param string opt: type of statistic: "N": number of elements. "mean": mean. "median": median. "min": minimum. "max": maximum. "stddev": standard deviation. [value]: percentile (a value between 0 and 100).

   Return value
      A scalar value.
  
   Description
      Calculates the relevant statistics for an input one-dimensional
      distribution. Missing values are ignored.
  
   Caveats
  
   Reference
  
   Modification history
      20140526-A_righ_ma: written.
  
.. function:: lognormal_dist(nn:numeric, dg:numeric, sig[1]:numeric, darr[*]:numeric)

   :param numeric nn: particle number concentration, can be a scalar or 1-D array.
   :param numeric dg: median diameter, same dimensionality of nn
   :param numeric sig: geometric standard deviation, a scalar
   :param numeric darr: array of diameters.

   Return value
      An array of type float, with the same dimensionality of nn, plus the darr
      dimension on the right, and with the same units of nn.
  
   Description
      Calculate a lognormal distribution given the three paramters and an array
      of diameters.
  
   Caveats
      dg and darr must have the same units.
  
   Reference
      Seinfeld and Pandis, Atmospheric chemistry and physics, JohnWiley & Sons,
      New York, US, 1998.
  
   Modification history
      20130528-A_righ_ma: written.
  
