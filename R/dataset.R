#' @title survival data
#' @docType data
#' @description A survival dataset related the long2 dataset, with different numeric and categorical covariate
#' @usage data(new_surv2)
#' @format a tibble of 13 columns and 1000 observations,
#' \describe{
#' \item{id}{id value for subjects}
#' \item{status}{survival status}
#' \item{time}{ survival time}
#' \item{visit}{visit time of longitudinal measurements}
#' \item{x1,x2,...,x7}{ different numeric and categorical variable}
#'
#' }
#'
"new_surv2"
#' @title longitudinal- survival dataset
#' @docType  data
#' @description A longitudinal dataset with single marker , with different numeric and categorical covariate
#' @usage data(new_long2)
#' @format a tibble of 13 columns and 5639 observations,
#' \describe{
#' \item{id}{id value for subjects}
#' \item{status}{survival status}
#' \item{time}{ survival time}
#' \item{y}{longitudinal marker}
#' \item{visit}{visit time of longitudinal measurements}
#' \item{x1,x2,...,x7}{ different numeric and categorical variable}
#'
#' }
#'
"new_long2"
