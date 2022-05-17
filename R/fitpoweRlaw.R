#' orbde
#'
#' This function returns the result of orbital decomposition. All of the
#' metrics associated with orbital decomposition are returned in a list.
#' @param <data_seq> character or integer time series
#' @keywords orbital decomposition
#' @export
#' @examples
#' orbde(guastello)
#'
orbde <- function(data_seq) {
  # Returns a table  with the final orbital decomposition analysis:
  #  - String length, C
  #  - Trace of the transition matrix for the string length C, M^C
  #  - Topological entropy, H_T (this is better than H_S for these purposes;
  #    use minimum of H_T)
  #  - Shannon entropy, H_S
  #  - D_L - dimensionality via Lyapunov
  #  - chi-square
  #  - degrees of freedom, df
  #  - p-value, p
  #  - Number of possible sub-strengths of length C, N* = length + 1 - C
  #  - phi-square (phi-sq = chi-sq/N*)
  # This all follows Guastello, Chapter 21, in Guastello & Gregson, 2011,
  #  _Nonlinear Dynamical Systems Analysis for the Behavioral Sciences Using
  #  Real Data_.
  #
  # 52 codes should be sufficient for any application; it is the limit on the
  #  ORBDE software as well.
  require(data.table)  # Necessary to allow parallelizing to work
  # data.table::shift() is used here.
  require(glue)

  # First, recast everything into single character codes:
  unique(data_seq) -> uni_seq
  length(uni_seq) -> n_codes
  if(n_codes > 52) {
    cat("Cannot do more than 52 codes\n") # See below.
    return(NULL)
  }
  c(LETTERS, letters)[1:n_codes] ->
    uni_rep                         # c(LETTERS, letters) has 52 elements
  uni_seq -> names(uni_rep)
  uni_rep[data_seq] -> data_seq

  # For now, remove all missing data from the sequence. In the case where
  #  the missing data is due to uncoded data, this will cause problems,
  #  but so be it.
  if(any(is.na(data_seq)) == TRUE) {
    data_seq[-which(is.na(data_seq))] -> data_seq
  }

  unique(data_seq) -> uni_seq
  length(uni_seq) -> n_codes

  table(data_seq) -> freqs
  # May want this for the previous line:
  # table(data_seq, useNA = "always") -> freqs

  glue_collapse(data_seq) -> coll_seq
  # Begin processing data: For C = 1
  1 -> C
  length(data_seq) -> seq_len -> N_star

  shift(data_seq, 1) -> shifted_seq
  length(which(data_seq == shifted_seq)) -> recurs

  if(recurs > 0) {
    unique(data_seq[data_seq == shifted_seq]) -> repeats




    repeats[!is.na(repeats)] -> repeats




    length(repeats) -> trMC
    if (trMC > 0) {
      log2(trMC) / C -> H_T
    } else {
      -Inf -> H_T
    }
    exp(H_T) -> D_L
  } else {
    0 -> trMC
    -Inf -> H_T
    0 -> D_L
  }
  n_codes - 1 -> dof # For the singlets, this is true;

  table(data_seq) -> freqs  # I think this is a repeated command,
  #  but everything breaks if I remove it.
  freqs / seq_len -> p_obs_keep  # This gets used for C > 1 calculations
  rep(seq_len / n_codes,
      n_codes) -> F_exp
  freqs / seq_len -> p_obs_tab
  - 1 * sum(p_obs_tab * log(p_obs_tab)) -> H_S

  sum(freqs * log(freqs / F_exp)) * 2 ->
    chi_sq

  dchisq(chi_sq, dof) -> p
  chi_sq / N_star -> phi_sq

  data.frame("C" = 1,
             "trM" = trMC,
             "Ht" = H_T,
             "Dl" = D_L,
             "chi^2" = chi_sq,
             "df" = dof,
             "N*" = N_star,
             "Phi^2" = phi_sq,
             "Hs" = H_S,
             "p" = p) -> OD_tab

  # Processing for 1 < C < N_star; while trace(C^M) > 0
  seq_len -> N_star
  for(len in 1:(length(data_seq) / 2)) {    # This is the ORBDE default
    # Identify all recurrences of length C
    C + 1 -> C
    N_star - 1 -> N_star

    # I suspect that this next part is the part that takes time
    vector(mode = "character",
           length = N_star) -> current_seq
    C -> end_index
    # The old way: routine is >10 times slower with this:
    #    for (index in 1:N_star) {
    #      paste0(data_seq[index:end_index],
    #             collapse = "") -> current_seq[index]
    #      end_index + 1 -> end_index
    #    }

    # New way: routine is much faster than with the old way
    for(index in 1:N_star) {
      substr(coll_seq, index, end_index) -> current_seq[index]
      end_index + 1 -> end_index
    }

    shift(current_seq, C) -> shifted_seq

    # Can we do a loop here and skip the which?
    # E.g.: if(current_seq == shifted_seq) {...}
    which(current_seq == shifted_seq) -> repeat_nums
    # I doubt it is faster.

    if (length(repeat_nums) > 0) {
      length(repeat_nums) -> recurs

      unique(current_seq[repeat_nums]) -> repeats  # Which codes are repeated?



      repeats[!is.na(repeats)] -> repeats # Just to be sure




      length(repeats) -> trMC                   # How many repeated codes
      #  are there?
      log2(trMC) / C -> H_T                  # Topological entropy
      exp(H_T) -> D_L                        # Lyapunov dimension

      # The number of unique repeated codes:
      table(current_seq) -> repeated_codes -> F_obs_tab
      length(repeated_codes[repeated_codes > 1]) -> dof
      # Shannon entropy. Must do this before collapsing codes:
      sum((F_obs_tab / N_star) *
            (log(N_star/F_obs_tab))) -> H_S

      F_obs_tab[F_obs_tab > 1] -> F_rep_obs_tab  # Repeated codes

      # Frequency expected
      length(F_rep_obs_tab) -> n_rep
      rep(1, n_rep) -> F_exp           # For all the repeats and

      # Is there a faster way here?
      for (i in 1:n_rep) {
        for (j in 1:C) {
          substr(names(F_rep_obs_tab)[i], j, j) -> fn
          p_obs_keep[fn] * F_exp[i] -> F_exp[i]
        }
      }
      F_exp * N_star -> F_exp

      # Guastello eq. 21.6:
      #  chi-squared = 2 * sum( F_obs * ln(F_obs/F_expected) )
      sum(F_rep_obs_tab * log(F_rep_obs_tab / F_exp)) * 2 -> chi_sq
      abs(N_star - sum(F_rep_obs_tab)) -> singles
      # Need to do something for all the non-recurrent sequences, as they
      #  contribute here too.
      if(singles > 0.01 ) {        # 0.01 is capricious, but probably good
        chi_sq + 2 *singles * log(singles / (N_star - sum(F_exp))) -> chi_sq
      }

      # Now for p-value and phi-squared (Guastello eq 21.7):
      dchisq(chi_sq, dof) -> p
      chi_sq / N_star -> phi_sq

      if(!is.na(trMC)) {
        if(trMC != 0) {
          data.frame("C" = C,
                     "trM" = trMC,
                     "Ht" = H_T,
                     "Dl" = D_L,
                     "chi^2" = chi_sq,
                     "df" = dof,
                     "N*" = N_star,
                     "Phi^2" = phi_sq,
                     "Hs" = H_S,
                     "p" = p) -> OD_line
          rbind(OD_tab, OD_line) -> OD_tab
        }
      }
    }
  }

  if(OD_tab$trM[1] == 0) {
    OD_tab[-1,] -> OD_tab
  }

  return(OD_tab)
}

#' IPL_fit
#'
#' This function returns the burstiness coefficient for each of the
#' unique codes in the time series.
#' @param <ts> character time series
#' @keywords burstiness
#' @export
#' @examples
#' IPL_fit(guastello, C)
#'
IPL_fit <- function(ts,       # data sequence
                    C = 1) {  # subsequences length
  # (from orbital decomposition, perhaps)
  # Create the subsequences

  unique(ts) -> num_codes

  matrix(data = 0,
         nrow = length(C),
         ncol = 3) -> fits_mat
  c("alpha", "R2", "H_N") -> colnames(fits_mat)

  # Need a list of codes
  code_list <- function(vec) {
    # This should never be true; all video has something in it
    if(is.na(vec[1]) == TRUE) {
      return(NULL)
    }
    if(length(which(is.na(vec[length(vec)]))) == 0) {
      return(vec[1:length(vec)])
    }
    return(vec[1:(min(which(is.na(vec)))-1)])
  }

  length(ts) -C + 1 -> len
  vector(mode = "character", length = len) ->
    sub_seq
  for (seq_start in 1:len) {
    paste0(ts[seq_start:(seq_start + C - 1)],
           collapse = "") ->
      sub_seq[seq_start]
  }

  table(sub_seq) -> freqs

  # Compute the entropy; see note following this chunk
  sum(freqs) -> Np1mC   # N (string length) plus 1 minus C
  -sum((freqs/Np1mC * log(freqs/Np1mC, base = 2))) /
    (Np1mC * num_codes^(-C) * log(num_codes^C, 2)) ->
    Entropy

  # Fit the data with a power law
  displ$new(unname(freqs)) -> m_pl
  if (length(unique(m_pl$dat)) > 2) {
    estimate_xmin(m_pl) -> est

      # The next occasionally happens; I think it is when the distribution
      #  is flat.
      if (m_pl$getXmin() > length(m_pl$internal$freq)) {
        1 -> est
      }
      m_pl$setXmin(est)

      if (!is.null(m_pl$getPars())) {
        # Get the exponent
        m_pl$getPars() -> Shape

        # Calculate R^2
        # Should this be used?
        cor(log(m_pl$internal$values), log(m_pl$internal$cum_n),
            method = "pearson"
        )^2 -> R2
      }
   }
  return(c(Shape, R2, Entropy))
}

