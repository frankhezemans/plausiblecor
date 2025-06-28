library(EMC2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(usethis)

design_forstmann <- EMC2::design(
  formula = list(
    v ~ lM, # mean drift coefficient: varying by S-R match
    A ~ 1,  # across-trial variability in start point: shared across conditions
    B ~ E,  # distance from start point to threshold: varying by speed emphasis
    t0 ~ 1, # non-decision time: shared across conditions
    sv ~ lM  # across-trial variability in drift coefficient: varying by S-R match
  ),
  model = EMC2::LBA,
  data = EMC2::forstmann,
  contrasts = list(
    # contrast matrices for effects of S-R match on v and sv
    v = list(
      lM = matrix(
        data = c(-1/2, 1/2),
        dimnames = list(NULL, "diff")
      )
    ),
    sv = list(
      lM = matrix(
        data = c(1/2, -1/2),
        dimnames = list(NULL, "diff")
      )
    ),
    # contrast matrix for effect of speed emphasis on distance to threshold:
    # First contrast represents the difference between speed and the other
    # conditions, matching the key test of Forstmann et al. (2008).
    # Second contrast represents the difference between accuracy and neutral
    # conditions, so that separate B values can be identified for each of the
    # three conditions.
    B = list(
      E = matrix(
        data = c(
          1/2, -1/4, -1/4,
          0, 1/2, -1/2
        ),
        ncol = 2,
        dimnames = list(NULL, c("speed_vs_avg", "neut_vs_acc"))
      )
    )
  ),
  # function to determine whether a given accumulator matches the stimulus
  matchfun = function(x) x$S == x$lR,
  # intercept of across-trial variability in drift coefficient fixed to satisfy
  # scaling constraint of LBA model
  constants = c(sv = log(1))
)

# Note the resulting design matrix for parameter B:
# $B
# E B B_Espeed_vs_avg B_Eneut_vs_acc
# speed 1            0.50            0.0
# neutral 1           -0.25            0.5
# accuracy 1           -0.25           -0.5
#
# Hence, the column "B" represents the intercept; the column "B_Espeed_vs_avg"
# represents the effect of speed relative to both neutral and accuracy; and the
# column "B_Eneut_vs_acc" represents the effect of neutral relative to accuracy.
# For "B_Espeed_vs_avg", more negative values correspond to a lower B parameter
# for speed relative to both neutral and accuracy. For "B_Eneut_vs_acc", more
# negative values correspond to a lower B parameter for neutral relative to
# accuracy.

prior_forstmann <- EMC2::prior(
  design = design_forstmann,
  type = "standard",
  pmean = c(
    v = 1.5, v_lMdiff = 2, sv_lMdiff = 0.5, A = log(0.5), t0 = log(0.15),
    B = log(0.5), B_Espeed_vs_avg = -0.5, B_Eneut_vs_acc = -0.25
  ),
  psd = c(
    v = 1.5, v_lMdiff = 2, sv_lMdiff = 0.5, A = 0.5, t0 = 0.25,
    B = 0.25, B_Espeed_vs_avg = 0.25, B_Eneut_vs_acc = 0.25
  )
)

emc_forstmann <- EMC2::make_emc(
  data = EMC2::forstmann,
  design = design_forstmann,
  type = "standard",
  n_chains = 4,
  prior_list = prior_forstmann
)

# NB following will take anywhere between a couple of minutes to several hours
# to run, depending on the capabilities of your machine
fit_forstmann <- EMC2::fit(
  emc = emc_forstmann,
  iter = 1500#,
  # fileName = "data-raw/Forstmann_EMC2_fit.RData"
)

# diagnostics for group-level parameters:
# R-hats all practically equal to 1; effective sample sizes well over 1000 per
# parameter.
summary(fit_forstmann, selection = "mu")
summary(fit_forstmann, selection = "sigma2")

# group-level hypothesis test: effect of speed condition (relative to both
# neutral and accuracy conditions) on distance to threshold.
# extremely strong evidence in favour of an effect: distance to threshold is
# strongly reduced in the speed condition.
hypothesis(
  emc = fit_forstmann,
  parameter = "B_Espeed_vs_avg",
  selection = "mu"
)

# extract posterior samples of participant-level parameters, and wrangle to
# obtain measures of speed effect on response caution, as in Forstmann et al.
Forstmann_LBA <- parameters(
  x = fit_forstmann,
  selection = "alpha"
) %>%
  # intermediate selection of columns
  dplyr::select(
    dplyr::all_of(c("subjects", "A")),
    dplyr::starts_with("B")
  ) %>%
  dplyr::mutate(
    subjects = factor(subjects)
  ) %>%
  # add sample ID column
  dplyr::mutate(
    .draw = rep(
      seq_len(4 * 1500),
      times = length(levels(EMC2::forstmann$subjects))
    ),
    .before = subjects
  ) %>%
  # transform parameter values
  dplyr::mutate(
    A = exp(A),
    B_speed = exp(B + (1/2) * B_Espeed_vs_avg),
    B_neutral = exp(B - (1/4) * B_Espeed_vs_avg + (1/2) * B_Eneut_vs_acc),
    B_accuracy = exp(B - (1/4) * B_Espeed_vs_avg - (1/2) * B_Eneut_vs_acc),
    .keep = "unused"
  ) %>%
  # pivot in long format by speed manipulation
  tidyr::pivot_longer(
    cols = dplyr::all_of(c(
      "B_speed", "B_neutral", "B_accuracy"
    )),
    names_to = "condition",
    names_prefix = "B_",
    values_to = "B"
  ) %>%
  # compute measure of response caution used in Forstmann et al. (2008):
  # ratio of LBA threshold (b) and starting point noise (A).
  # we have estimates of A and the *distance from* A to b, termed B:
  # B = b - A. hence, threshold b = A + B.
  dplyr::mutate(
    threshold = A + B,
    caution = threshold / A,
    .keep = "unused"
  ) %>%
  dplyr::select(
    -dplyr::all_of("threshold")
  ) %>%
  # spread back into wide format
  tidyr::pivot_wider(
    values_from = "caution",
    names_from = "condition"
  ) %>%
  # lastly, compute the desired measure of the effect of the speed condition on
  # response caution
  dplyr::mutate(
    caution_effect_speed = speed - ((neutral + accuracy) / 2),
    .keep = "unused"
  )

usethis::use_data(Forstmann_LBA, overwrite = TRUE, compress = "xz")
