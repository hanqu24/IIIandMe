#' Generating all possible haplotypes & positions of crossovers
#'
#' @param data preprocessing data
#' @param npol pollen numbers
#' @param filterGenoError filtering genotyping errors
#'
#' @return lists of possible haplotypes & crossovers
#' @export
#'
#' @examples
#' library(IIIandme)
#' input<- PreProcessing(sample)
#' HapCo(input, 5)

HapCo <- function(data, npol, filterGenoError = F) {
  library(gtools)
  choice <- gtools::combinations(npol, 3, 5:(npol + 4))
  res <- hap3fromLot(data, choice, filter = filterGenoError)
  return(res)
}

hap3fromLot <- function(cpn, choice, filter) {
  library(Hapi)
  hap <- list()
  co <- list()
  for (a in 1:nrow(choice)) {
    print(nrow(cpn))
    ### Filter NAs
    pn <- basicFilterFun(cpn[, choice[a, ]], 3)
    print(nrow(pn))

    if (filter == TRUE) {
      pn <- hapiFilterError(pn)
    }

    print("pollens are ------")
    print(colnames(pn))
    print("--------")
    p1 <- as.numeric(pn[, 1])
    p2 <- as.numeric(pn[, 2])
    p3 <- as.numeric(pn[, 3])

    inferredParent <- phase3pollens(p1, p2, p3)
    hap[[a]] <- inferredParent[[1]]
    rownames(hap[[a]]) <- rownames(pn)
    if (sum(is.na(inferredParent[[2]][, 3])) == nrow(inferredParent[[2]])) {
      co[[a]] <- data.frame(inferredParent[[2]], pos = rep(NA, nrow(inferredParent[[2]])))
    } else {
      co[[a]] <- data.frame(inferredParent[[2]], pos = rownames(pn)[inferredParent[[2]][, 3]])
    }
    co[[a]][, 1][which(co[[a]][, 1] == "pn1")] <- colnames(pn)[1]
    co[[a]][, 1][which(co[[a]][, 1] == "pn2")] <- colnames(pn)[2]
    co[[a]][, 1][which(co[[a]][, 1] == "pn3")] <- colnames(pn)[3]

    out <- list(hap, co)
  }
  return(out)
}



phase3pollens <- function(pn1, pn2, pn3) {
  ### two complete chromatids with exactly same genotype ###
  if (prod(pn1 == pn2) == 1) {
    print("Pollen #1 & #2 are same")
    hap_1 <- pn1
    hap_2 <- flipFun(pn1)

    hap_loci <- list(cbind(hap_1, hap_2), FindThirdCrossover(pn3, hap_1, "pn3", "pn1", "pn2"))
    rownames(hap_loci[[1]]) <- names(pn1)
    return(hap_loci)
  }
  if (prod(pn1 == pn3) == 1) {
    print("Pollen #1 & #3 are same")
    hap_1 <- pn1
    hap_2 <- flipFun(pn1)

    hap_loci <- list(cbind(hap_1, hap_2), FindThirdCrossover(pn2, hap_1, "pn2", "pn1", "pn3"))
    rownames(hap_loci[[1]]) <- names(pn1)
    return(hap_loci)
  }
  if (prod(pn2 == pn3) == 1) {
    print("Pollen #2 & #3 are same")
    hap_1 <- pn2
    hap_2 <- flipFun(pn2)

    hap_loci <- list(cbind(hap_1, hap_2), FindThirdCrossover(pn1, hap_1, "pn1", "pn2", "pn3"))
    rownames(hap_loci[[1]]) <- names(pn1)
    return(hap_loci)
  }

  ### two complete chromatids with exactly complementary genotype ###
  if (prod(pn1 == flipFun(pn2)) == 1) {
    print("Pollen #1 & fliped #2 are same")
    hap_1 <- pn1
    hap_2 <- pn2

    hap_loci <- list(cbind(hap_1, hap_2), FindThirdCrossover(pn3, hap_1, "pn3", "pn1", "pn2"))
    rownames(hap_loci[[1]]) <- names(pn1)
    return(hap_loci)
  }
  if (prod(pn1 == flipFun(pn3)) == 1) {
    print("Pollen #1 & fliped #3 are same")
    hap_1 <- pn1
    hap_2 <- pn3

    hap_loci <- list(cbind(hap_1, hap_2), FindThirdCrossover(pn2, hap_1, "pn2", "pn1", "pn3"))
    rownames(hap_loci[[1]]) <- names(pn1)
    return(hap_loci)
  }
  if (prod(pn2 == flipFun(pn3)) == 1) {
    print("Pollen #2 & fliped #3 are same")
    hap_1 <- pn2
    hap_2 <- pn3

    hap_loci <- list(cbind(hap_1, hap_2), FindThirdCrossover(pn1, hap_1, "pn1", "pn2", "pn3"))
    rownames(hap_loci[[1]]) <- names(pn1)
    return(hap_loci)
  }

  # print("No complete chromatids have been sampled")

  ### Majority voting walking ###

  ### Check if there is common region(s) across pollens

  pn <- rbind(pn1, pn2, pn3)
  isCommon <- findCommonFun(pn)

  # Need a common region that is at least 3 markers long
  if (length(isCommon) == 0) {
    ind <- 0
    ii <- 1
    while (ii <= 3) {
      pn[ii, ] <- flipFun(pn[ii, ]) # Flip one pollen
      isCommon <- findCommonFun(pn)
      if (length(isCommon) == 0) {
        pn[ii, ] <- flipFun(pn[ii, ]) # Does not create common region; flip back this pollen
        ii <- ii + 1
      } else {
        ind <- 1
        message(sprintf("Common region found after fliping %s\n", ii))
        break
      }
    }

    if (ind == 1) {
      isCommon <- findCommonFun(pn)
    } else {
      stop("Cannot find common region even after fliping")
    }
  }


  ### Start walking from the common region (backbone)

  start <- isCommon[1]
  end <- isCommon[2]

  message(sprintf("Common region found between: %s\n", isCommon))

  backbone <- matrix("b", 3, (end - start + 1))

  # print("Initial backbone set up:")
  # print(backbone)

  if (end != ncol(pn)) {
    ### Walking towards right
    # print("Walking towards right.")
    v0 <- pn[, end]
    ch0 <- c("b", "b", "b")
    for (i in (end + 1):ncol(pn)) {
      ch.temp <- walking3Fun(v0, ch0, pn[, i])
      backbone <- cbind(backbone, ch.temp)
      v0 <- pn[, i]
      ch0 <- ch.temp
    }
    print("..................................................................................")
    print("Walking towards right done.")
    print("..................................................................................")
  }

  if (start != 1) {
    ### Walking towards left
    # print("Walking towards left.")
    v0 <- pn[, start]
    ch0 <- c("b", "b", "b")
    for (i in rev(1:(start - 1))) {
      ch.temp <- walking3Fun(v0, ch0, pn[, i])
      backbone <- cbind(ch.temp, backbone)
      v0 <- pn[, i]
      ch0 <- ch.temp
    }
    print("..................................................................................")
    print("Walking towards left done.")
    print("..................................................................................")
  }

  p1l <- crossoverCountFun(backbone[1, ])
  p2l <- crossoverCountFun(backbone[2, ])
  p3l <- crossoverCountFun(backbone[3, ])

  message(sprintf("The number & position of crossovers on pollen #1 is: %s\n", p1l))
  message(sprintf("The number & position of crossovers on pollen #2 is: %s\n", p2l))
  message(sprintf("The number & position of crossovers on pollen #3 is: %s\n", p3l))

  l1 <- cbind(data.frame(pollen = rep("pn1", nrow(p1l))), p1l)
  l2 <- cbind(data.frame(pollen = rep("pn2", nrow(p2l))), p2l)
  l3 <- cbind(data.frame(pollen = rep("pn3", nrow(p3l))), p3l)


  if (which.min(c(p1l[1, 1], p2l[1, 1], p3l[1, 1])) == 1) {
    hap_1 <- ifelse(backbone[1, ] == "b", pn1, flipFun(pn1))
    hap_2 <- flipFun(hap_1)
  } else if (which.min(c(p1l[1, 1], p2l[1, 1], p3l[1, 1])) == 2) {
    hap_1 <- ifelse(backbone[2, ] == "b", pn2, flipFun(pn2))
    hap_2 <- flipFun(hap_1)
  } else {
    hap_1 <- ifelse(backbone[3, ] == "b", pn3, flipFun(pn3))
    hap_2 <- flipFun(hap_1)
  }

  hap_loci <- list(cbind(hap_1, hap_2), rbind(l1, l2, l3))
  rownames(hap_loci[[1]]) <- names(pn1)
  return(hap_loci)
}
