#' Interpret Generalised Additive Model Output
#'
#' @param mod model output of class `gam`
#' @param words integer, indicating different contrasting words. 0 = the numerical
#'  output; 1 = 'more' or 'less'; 2 = 'greater' or 'less'; 3 = 'higher' or 'lower'
#' @param rounding integer, indicating the number of decimals
#' @param coeff character, specifying the coefficient name as it appears in the
#'  model output. Defaults to "(Intercept)".
#' @param sig logical, if TRUE return the statistical significance associated
#'  with the `coeff`. words = 1 returns "not significant" or "significant".
#'  `words = 0` returns the numeric output.
#'
#' @returns character or numeric depending on `words`
#' @export
#'
#' @examples
#' mod1 <- lm(speed ~ dist, data = cars)
#' interpret
interpret_gam <- function(mod,
                          words = 1,
                          coeff = "(Intercept)",
                          sig = FALSE,
                          rounding = 2) {

  words <- as.character(words)

  if(sig){
    out <- switch(words,
                  "0" = round(unname(summary(mod)$p.table[,4][coeff]),rounding),
                  "1" = ifelse(summary(mod)$p.table[,4][coeff] >0.05,
                             "not significant",
                             "significant"))
    }else{

  out <- switch(words,
                "0" = round(unname(mod$coefficients[coeff]),rounding),
                "1" = ifelse(mod$coefficients[coeff] > 0, "more", "less"),
                "2" = ifelse(mod$coefficients[coeff] > 0, "greater", "fewer"),
                "3" = ifelse(mod$coefficients[coeff] > 0, "higher", "lower"),
                "4" = ifelse(mod$coefficients[coeff] > 0, "after", "before"),
                "5" = ifelse(mod$coefficients[coeff] > 0, "later", "earlier"))
  }

  return(unname(out))

}
