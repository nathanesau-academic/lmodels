#' @title Ogive
#' @description Computes the ogive value for a given x
#' @param x The value to pass to the ogive function
#' @param cdf A cdf function to use (use ecdf on the data)
#' @param breaks The breaks for the ogive function
#' @return The ogive value (between 0 and 1) for a given x and cdf function.
#' @examples x <- c(8,9,13,15,16,22,23,27,29,35)
#' c <- c(0,10,15,20,30,40)
#' f <- ecdf(x)
#' plot(function(y) ogive(y,f,c), 0, 50, main="Ogive", ylab="F(x)")
#' @details The cdf is computed before the ogive function is called so
#'          it doesn't have to be recomputed every time ogive(...) is called.
#' @export
ogive <- function(x, cdf, breaks) {
  f <- function(x) {
    if(x<breaks[1])        return(0)
    if(x>=tail(breaks,1))  return(1)

    cindex <- tail(which(breaks<=x),1)
    c0 <- breaks[cindex]
    c1 <- breaks[cindex+1]

    cdf(c0) + (cdf(c1)-cdf(c0))/(c1 - c0) * (x - c0)
  }
  sapply(x,f)
}

#' @title Risk Set
#' @description Computes the risk set for a given vector x
#' @param x The data
#' @param d The left truncation points (optional)
#' @param u The right censoring points (optional)
#' @return A data frame containing the unique values of x (y),
#'         the counts for each unique value (s), and the risk
#'         set (r)
#' @examples x <- c(1.0,1.3,1.5,1.5,2.1,2.1,2.1,2.8)
#' rs <- RiskSet(x)
#' @details If d and u are provided, they should be the same lenggth as
#' x. Entries of d, u and x can be set to NA where applicable.
#' @export
RiskSet <- function(x,d=NA,u=NA) {
  y <- as.numeric(names(table(x)))
  s <- as.numeric(table(x))
  r <- sapply(y, function(m) sum(x[!is.na(x)]>=m) +
                sum(u[!is.na(u)]>=m) -
                sum(d[!is.na(d)]>=m))
  data.frame(y=y,s=s,r=r)
}

#' @title Kaplan Meier Product Limit Estimator
#' @description Survival function estimate using Kaplan Meier method
#' @param x The value to pass to the Kaplan Meier function
#' @param rs is the output from RiskSet(...)
#' @return The estimate of the survival function at a given value of x
#' @examples x <- c(1.0,1.3,1.5,1.5,2.1,2.1,2.1,2.8)
#' rs <- RiskSet(x)
#' plot.stepfun(ecdf(x), verticals=TRUE, pch='',
#'              ylab="F(x)", main="Empirical CDF")
#' plot(function(y) NelsonAalen(y,rs), 0, 3,
#'      type='s', ylab="H(x)", main="Nelson-Aalen estimate of H(x)")
#' @export
KaplanMeier <- function(x,rs) {
  f <- function(x) {
    if(x<rs$y[1])   return(1)
    yindex <- tail(which(rs$y<=x),1)
    prod((rs$r[1:yindex] - rs$s[1:yindex])/rs$r[1:yindex])
  }

  sapply(x,f)
}

#' @title Nelson Aalen Estimator
#' @description Cumulative hazard function estimate using Nelson Aalen method
#' @param x The value to pass to the Nelson Aalen function
#' @param rs is the output from RiskSet(...)
#' @param survival If TRUE return estimate of survival function, otherwise return estimate of cumulative hazard function
#' @return The estimate  of the cumulative hazard function at a given value of x
#' @details To get estimate of survival function take exp(-H(x)) where H(x) is the cumulative hazard rate
#' @examples x <- c(1.0,1.3,1.5,1.5,2.1,2.1,2.1,2.8)
#' rs <- RiskSet(x)
#' plot(function(y) NelsonAalen(y,rs), 0, 3,
#'      type='s', ylab="H(x)", main="Nelson-Aalen estimate of H(x)")
#' @export
NelsonAalen <- function(x,rs,survival=FALSE) {
  f <- function(x) {
    if(x<rs$y[1]) {
      if(survival)  return(1)
      else          return(0)
    }
    yindex <- tail(which(rs$y<=x),1)
    H <- sum(rs$s[1:yindex]/rs$r[1:yindex])
    ifelse(survival, exp(-H), H)
  }

  sapply(x,f)
}
