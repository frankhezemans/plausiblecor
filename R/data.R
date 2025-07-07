#' Forstmann et al. (2008) fMRI contrasts
#'
#' Estimated fMRI contrasts from Forstmann et al. (2008), representing the
#' percentage of BOLD signal change in two brain regions: the right anterior
#' striatum (caudate nucleus) and the right pre-supplementary motor area
#' (pre-SMA).
#'
#' The data were collected while participants performed a speeded perceptual
#' decision-making task. On each trial, participants judged whether a cloud of
#' moving dots appeared to move to the left or to the right. Before each trial,
#' an **instructional cue** indicated the emphasis for that trial: either to
#' respond as **quickly** as possible (speed condition), as **accurately** as
#' possible (accuracy condition), or to adopt a **balanced** speed–accuracy
#' trade-off (neutral condition). These cues were assumed to modulate response
#' caution on a trial-by-trial basis.
#'
#' The fMRI contrasts reflect results from a group-level conjunction analysis
#' that identified voxels more active in both the *speed > accuracy* and *speed > neutral*
#' contrasts, using a cluster-level correction for multiple comparisons at
#' *P* < 0.001. Brain regions of interest were thus **functionally defined** based on
#' this contrast-conjunction, and included contiguous voxel clusters in the
#' right anterior striatum and right pre-SMA that exceeded the critical z-threshold.
#'
#' To compute the percentage signal change, the time course of the BOLD response
#' was extracted from the activated voxel clusters in each region for each participant.
#' The response was measured relative to a baseline defined by null events and
#' averaged across the 4–8 second window following presentation of the **instructional cue**.
#' These values reflect increased neural activation in the speed condition relative
#' to both the neutral and accuracy conditions.
#'
#' @format ## `Forstmann_fMRI`
#' A data frame with 19 rows and 3 columns:
#' \describe{
#'   \item{subjects}{Factor variable: de-identified participant ID codes.}
#'   \item{striatum}{Numeric variable: percentage BOLD signal change in the right caudate nucleus.}
#'   \item{pre_sma}{Numeric variable: percentage BOLD signal change in the right pre-SMA.}
#' }
#'
#' @source <https://osf.io/download/f7553/>
#'
#' @references
#' Forstmann, B. U., Dutilh, G., Brown, S., Neumann, J., von Cramon, D. Y.,
#' Ridderinkhof, K. R., & Wagenmakers, E.-J. (2008). Striatum and pre-SMA
#' facilitate decision-making under time pressure. *Proceedings of the National
#' Academy of Sciences*, 105(45), 17538–17542. \doi{10.1073/pnas.0805903105}
"Forstmann_fMRI"

#' Forstmann et al. (2008) response caution contrast
#'
#' Estimated effect on the response caution parameter from a linear ballistic
#' accumulator (LBA) model, fit to data from Forstmann et al. (2008).
#'
#' Participants performed a speeded perceptual decision-making task in which they judged
#' whether a cloud of moving dots appeared to move left or right. Before each trial, an
#' instructional cue indicated the required emphasis: to respond as **quickly** as
#' possible (speed condition), as **accurately** as possible (accuracy condition),
#' or with a **balanced** speed–accuracy trade-off (neutral condition).
#' These cues were assumed to modulate response caution on a trial-by-trial basis.
#'
#' The data were fit using an LBA model (Brown & Heathcote, 2008) with a parameterisation
#' conceptually similar to that used by Forstmann et al. (2008). The LBA assumes evidence
#' accumulates linearly at a rate drawn from a normal distribution (mean `v`, SD `sv`)
#' until reaching a threshold `b`, at which point a response is made. The accumulation
#' process starts from a point drawn from a uniform distribution between 0 and `A`.
#' A non-decision time `t0` is added to account for perceptual and motor latencies.
#' The threshold was indirectly estimated via the non-negative quantity `B = b - A`.
#'
#' To test the hypothesis that the speed–accuracy manipulation influenced response
#' caution, `B` was allowed to vary across the three conditions (speed, neutral,
#' accuracy). To account for differences in response accuracy, the mean and
#' standard deviation of drift rate (`v` and `sv`) varied by response correctness.
#' These effects were modeled using sum-to-zero contrasts, where parameter values
#' that differ by levels of a categorical variable
#' (perceptual cue for `B`; response correctness for `v` and `sv`) are coded as
#' deviations from a grand mean (intercept). To satisfy the scaling constraint of
#' the LBA model (Donkin et al., 2009), the intercept of `sv` was arbitrarily fixed
#' to 1. The parameters `A` and `t0` were assumed to be shared across conditions.
#'
#' Model parameters were estimated using hierarchical Bayesian methods via the
#' EMC2 package (Stevenson et al., 2025), drawing 6000 posterior samples
#' (4 chains × 1500 iterations). Following Forstmann et al. (2008), response caution was
#' quantified as the ratio `b / A`. The reported contrast is the difference in response
#' caution between the speed condition and the mean of the neutral and accuracy
#' conditions, providing an estimate of the effect of speed emphasis on response caution.
#'
#' While the modelling approach here is conceptually similar to that of Forstmann et al.
#' (2008), the estimation method differs substantially (hierarchical Bayesian vs.
#' non-hierarchical maximum likelihood). Parameter estimates are therefore not expected
#' to match exactly.
#'
#' Note also that subsequent re-analyses suggest more complex model specifications may
#' better account for these data (see, e.g., the
#' [EMC2 tutorial chapter on the LBA model](https://bookdown.org/reilly_innes/EMC_bookdown/using-emc2-for-eams-ii---race-models.html)).
#'
#' @format ## `Forstmann_LBA`
#' A data frame with 114,000 rows and 3 columns:
#' \describe{
#'   \item{.draw}{Integer. Posterior sample number.}
#'   \item{subjects}{Factor. De-identified participant ID.}
#'   \item{caution_effect_speed}{Numeric. Estimated difference in response caution
#'   between the speed condition and the average of the neutral and accuracy conditions.}
#' }
#'
#' @source [EMC2::forstmann]
#'
#' @references
#' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice
#' response time: Linear ballistic accumulation. *Cognitive Psychology*, 57(3),
#' 153–178. \doi{10.1016/j.cogpsych.2007.12.002}
#'
#' Donkin, C., Brown, S. D., & Heathcote, A. (2009). The overconstraint of
#' response time models: Rethinking the scaling problem.
#' *Psychonomic Bulletin & Review*, *16*, 1129-1135. \doi{10.3758/PBR.16.6.1129}
#'
#' Forstmann, B. U., Dutilh, G., Brown, S., Neumann, J., von Cramon, D. Y.,
#' Ridderinkhof, K. R., & Wagenmakers, E.-J. (2008). Striatum and pre-SMA
#' facilitate decision-making under time pressure. *Proceedings of the National
#' Academy of Sciences*, 105(45), 17538–17542. \doi{10.1073/pnas.0805903105}
#'
#' Stevenson, N., Donzallaz, M. C., & Heathcote, A. (2025). *EMC2: Bayesian
#' Hierarchical Analysis of Cognitive Models of Choice*. R package version 3.1.1.
"Forstmann_LBA"
