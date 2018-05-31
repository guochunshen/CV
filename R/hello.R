#' Calculate coefficient of variation (CV)
#'
#' @param traits a vector of sampled trait values, minimum sample size is 10.
#'
#' @details
#' Given a vector of sampled values, this function will return a vector of CV by six estimators.
#'
#' @return
#' a vector of CV estimates
#'
#'
#' @references
#' Yang et al. (2018) How to accurately estimate intraspecific trait variation? Methods in Ecology and Evolution (submitted)
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

