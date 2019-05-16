HSSPScore <- function(PID,
                      MatchLength) {
  if (MatchLength <= 11L) {
    Score <- PID - 100L
  } else if (MatchLength > 11L && MatchLength <= 450L) {
    L <- -0.32 * (1 + exp(-MatchLength/1000))
    Score <- PID - (480L * MatchLength^L)
  } else {
    Score <- PID - 19.5
  }
  return(Score)
}