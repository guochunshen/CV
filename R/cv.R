#' Calculate coefficient of variation (CV)
#'
#' \code{cv} estimate CV by various CV estimators from a given set of samples
#'
#' @param traits a vector of sampled trait values, minimum sample size is 10.
#'
#' @details
#' Given a vector of sampled values, this function will return a vector of CV by six estimators.
#'
#' \describe{
#'  \item{CV1}{the most commonly used CV estimator, \eqn{CV1=sd_sample/mu_sample},
#'             where sd_sampe and mu_sample are the standard deviation and mean of the sample values. }
#'  \item{CV2}{bias corrected CV estimator under the assumption of normal trait distribution,
#'              \eqn{CV2=sd_sample/mu_sample(1+1/4/N), where N is the sample size.}}
#'  \item{CV3}{approximate CV estimator by Breunig (2001) without assumption of sample value distribution}
#'  \item{CV4}{approximate CV estimator by Bao (2009) without assumption of sample value distribution}
#'  \item{CV5}{composite CV esetiamtor, it is the arithmetic mean of CV3 and CV4, see details in Yang (2018)}
#'  \item{CV6}{composite CV esetiamtor, it is the arithmetic mean of CV2 and CV4, see details in Yang (2018)}
#' }
#'
#' @return
#' a vector of CV estimates
#'
#'
#' @references
#' \enumerate{
#'  \item Breunig, R. (2001) An almost unbiased estimator of the coefficient of variation. Economics Letters, 70 (1), 15–19.
#'  \item Bao, Y. (2009) Finite-sample moments of the coefficient of variation. Econometric Theory, 25 (01), 291–297.
#'  \item Yang et al. (2018) How to accurately estimate intraspecific trait variation? Methods in Ecology and Evolution (submitted)
#' }
#'
#'
#' @examples
#'
#' #simulated trait values
#' traits=rnorm(1000, 10, 1)
#'
#' #calculate the CV
#' cv(traits)
#'
cv=function(traits){

  if(any(is.na(traits))){
    warning("Missing values detected in samples and will be removed in the analysis.")
    traits=traits[!is.na(traits)]
  }

  N=length(traits)

  if(N<10)
    stop("sample size is too small, minimum sample size is 10.")

  #calcualte CV^2
  y_bar=mean(traits)
  s2_hat=var(traits)
  cv_2=s2_hat/y_bar^2
  cv_1=sqrt(cv_2)
  gamma_1=sum(((traits-y_bar)/s2_hat^0.5)^3)/N
  gamma_2=sum(((traits-y_bar)/s2_hat^0.5)^4)/N
  bias=cv_2^(3/2)/N*(3*cv_2^0.5-2*gamma_1)
  bias2=cv_1^3/N-cv_1/4/N-cv_1^2*gamma_1/2/N-cv_1*gamma_2/8/N
  cv1=sd(traits)/mean(traits)
  cv2=(1+1/4/N)*cv1
  cv3=sqrt(cv_2-bias)
  cv4=cv_1-bias2
  cv5=mean(cv3,cv4)
  cv6=mean(cv2,cv4)
  re=c(cv1,cv2,cv3,cv4,cv5,cv6)
  names(re)=paste("CV",1:6,sep="")
  return(re)
}

