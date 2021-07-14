#' compute.threshold.YI.AROC.sp function modified to retrieve sensitivity and specificity
#'
#' @param data https://rdrr.io/cran/ROCnReg/src/R/compute.threshold.YI.AROC.sp.R
#'
#' @return compute.threshold.YI.AROC.sp with sensitivity and specificity
#' @export
#'
#' @examples
aroc_compute<-function (object, newdata, parallel = c("no", "multicore", "snow"),
                 ncpus = 1, cl = NULL)
{
  if (class(object)[1] != "AROC.sp") {
    stop(paste0("This function can not be used for this object class: ",
                class(object)[1]))
  }
  names.cov <- all.vars(object$formula)[-1]
  if (!missing(newdata) && !inherits(newdata, "data.frame"))
    stop("Newdata must be a data frame")
  if (!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov,
                                                                     names(newdata)))))
    stop("Not all needed variables are supplied in newdata")
  if (missing(newdata)) {
    newdata <- cROCData(object$data, names.cov, object$group)
  }
  else {
    newdata <- na.omit(newdata[, names.cov, drop = FALSE])
  }
  p <- seq(0, 1, length = 500)
  np <- length(p)
  npred <- nrow(newdata)
  sigma0 <- summary(object$fit)$sigma
  data.d <- (object$data[object$data[, object$group] != object$tag.h,
  ])[!object$missing.ind$d, ]
  n1 <- nrow(data.d)
  pre.placement.values <- (data.d[, object$marker] - predict(object$fit,
                                                             newdata = data.d))/sigma0
  if (object$est.cdf.h == "normal") {
    u1 <- 1 - pnorm(pre.placement.values)
  }
  else {
    res0p <- object$fit$residuals/sigma0
    F0res <- ecdf(res0p)
    u1 <- 1 - F0res(pre.placement.values)
  }
  AROC <- numeric(np)
  for (i in 1:np) {
    AROC[i] <- sum(u1 <= p[i])/n1
  }
  difbb <- AROC - p
  FPF <- mean(p[which(difbb == max(difbb))])
  YI <- max(difbb)
  pred0 <- predict(object$fit, newdata = newdata)
  if (object$est.cdf == "normal") {
    csf0_inv <- qnorm(1 - FPF)
  }
  else {
    csf0_inv <- quantile(res0p, 1 - FPF, type = 1)
  }
  thresholds <- pred0 + sigma0 * csf0_inv
  res <- list()
  res$call <- match.call()
  res$newdata <- newdata
  res$thresholds <- thresholds
  res$YI <- YI
  res$FPF <- FPF
  res$specificity <- (1-p)
  res$sensitivity <- (AROC)
  res
}
