# ashr v2.0

## Major Changes

There are a few major changes in output and input that will likely
break existing dependencies. Here are the highlights

- the main output (lfsr, lfdr, etc) is rearranged into a dataframe,
  called `result`.

    - so, for example, the lfsr is now `a$result$lfsr` instead of
      `a$lfsr`

    - or, better, use the accessor function `get_lfsr(a)` to extract
      the lfsr etc

- I added accessor functions `get_lfsr`, `get_lfdr`, `get_pm`,
  `get_psd` etc to access the lfsr, lfdr, posterior mean and posterior
  standard deviation. Using these functions to access results will
  help ensure your code remains valid if I happen to change the
  internal structure of the results again (although not
  anticipated...)

- output `fitted.g` is renamed `fitted_g`, and `flash.data` becomes
  `flash_data` to make the whole package convention more
  consistent. Also `fit` becomes `fit_details`.

- function prefixes `comppost` and `compdens` replaced with
  `comp_post` and `comp_dens`, again for consistency.

- `nonzeromode` option is replaced with the option `mode` to specify
  mode. Or use `mode="estimate"` to estimate the mode.

- more flexible control of output. For example, you can say you want
  only the logLR output by specifying `outputlevel =
  c("lfsr","logLR")`, or only posterior mean by `outputlevel =
  c("PosteriorMean")`.
  
