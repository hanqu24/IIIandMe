#' Preprocessing & Basic functions
#'
#' @param data raw data as the sample data
#'
#' @return output of preprocessing data
#' @export
#'
#' @examples
#' library(IIIandMe)
#' PreProcessing(sample)

PreProcessing <- function(data) {
  library(tidyr)
  homo <- data %>% tidyr::separate(geno, c("A", "B", "C"), "")
  heter <- homo[which(homo$B != homo$C), ]
  heter <- heter[which(heter$ref != heter$alt), ]

  result <- cbind(heter[, 1:4], apply(heter[, 8:ncol(heter)], 2, function(x) ifelse(x == heter$ref, 1, 0)))
  return(result)
}

flipFun <- function(v) {
  v2 <- ifelse(v == 0, 1, 0) ### Complementary
  return(v2)
}

### Keep SNPs observed in at least n gamates ###
basicFilterFun <- function(v, n = 3) {
  v <- v[which(apply(v, 1, function(y) sum(!is.na(y))) >= n), ]
  return(v)
}

### Find the common region across pollens ###
findCommonFun <- function(v) {
  ll <- as.numeric(apply(v, 2, identicalFun))
  return(findLengthOfFirst1Fun(ll))
}

findAllCommonFun <- function(v) {
  ll <- as.numeric(apply(v, 2, identicalFun))
  return(crossoverCountFun(ll))
}


FindThirdCrossover <- function(pndiff, hap, ndiff, n1, n2) {
  c <- findAllCommonFun(t(cbind(pndiff, hap)))
  f <- cbind(pollen = rep(ndiff, nrow(c)), c)
  na <- data.frame(pollen = c(n1, n2), number = c(0, 0), locus = c(NA, NA))
  return(rbind(f, na))
}



### Check if all elements in vector are the same ###
identicalFun <- function(v) {
  return(length(unique(v)) == 1)
}


### Find the first string of 1's that has at 3 markers long ###
findLengthOfFirst1Fun <- function(v) {
  com <- NULL
  if (length(which(v == 1)) == 0) {
    print("Cannot find common regions")
    return(com)
  } else {
    start <- which(v == 1)[[1]]
    while (length(start) != 0) {
      if (start == 1) {
        if (length(which(v == 0)) != 0) {
          end <- which(v == 0)[[1]] - 1
        } else {
          end <- length(v)
        }
        if (end >= 3) {
          com <- c(start, end)
          return(com)
        } else {
          v[1:end] <- 9 ### mask the past region
          if (length(which(v == 1)) == 0) {
            print("Cannot find common region of length 3+")
            return(NULL)
          } else {
            start <- which(v == 1)[[1]]
          }
        }
      } else {
        v[1:(start - 1)] <- 9 ### mask the past region
        if (length(which(v == 0)) != 0) {
          end <- which(v == 0)[[1]] - 1
        } else {
          end <- length(v)
        }
        if ((end - start + 1) >= 3) {
          com <- c(start, end)
          return(com)
        } else {
          v[1:end] <- 9
          if (length(which(v == 1)) == 0) {
            print("Cannot find common region of length 3+")
            return(NULL)
          } else {
            start <- which(v == 1)[[1]]
          }
        }
      }
    }
  }

  print("Searching common regions done")
  return(com)
}


### Check pattern of genotypes at a location across 3 pollens ###
pattern3Fun <- function(v) {
  if (v[1] == v[2] & v[1] == v[3]) {
    return(1)
  } else if (v[1] == v[2] & v[1] != v[3]) {
    return(2)
  } else if (v[1] != v[2] & v[1] == v[3]) {
    return(3)
  } else if (v[1] != v[2] & v[1] != v[3] & v[2] == v[3]) {
    return(4)
  } else {
    return(0)
  }
}

flipChainFun <- function(v) {
  if (v != "b" & v != "c") {
    stop("Error")
  } else if (v == "b") {
    return("c")
  } else {
    return("b")
  }
}

crossoverCountFun <- function(v) {
  v <- as.numeric(as.factor(v))
  v1 <- v[-1]
  v2 <- v[-length(v)]
  vd <- v1 - v2
  if (prod(vd == 0) == 1) {
    return(data.frame(number = 0, locus = NA))
  } else {
    return(data.frame(number = sum(vd != 0), locus = which(vd != 0) + 1))
  }
}

### Walking to extend backbone ###
walking3Fun <- function(v1, ch1, v2) {
  pat.1 <- pattern3Fun(v1)
  pat.2 <- pattern3Fun(v2)
  if (prod(pat.1 == pat.2) == 1) {
    ch2 <- ch1
  } else {
    ch2 <- ch1
    pat <- sort(c(pat.1, pat.2))
    if (prod(pat == c(1, 2)) == 1) {
      ch2[3] <- flipChainFun(ch2[3])
    } else if (prod(pat == c(1, 3)) == 1) {
      ch2[2] <- flipChainFun(ch2[2])
    } else if (prod(pat == c(1, 4)) == 1) {
      ch2[1] <- flipChainFun(ch2[1])
    } else if (prod(pat == c(2, 3)) == 1) {
      ch2[1] <- flipChainFun(ch2[1])
    } else if (prod(pat == c(2, 4)) == 1) {
      ch2[2] <- flipChainFun(ch2[2])
    } else {
      ch2[3] <- flipChainFun(ch2[3])
    }
  }
  return(ch2)
}
