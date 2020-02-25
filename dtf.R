###############################################################################################
###############################################################################################
## Analyses of differential test functioning following Chalmers et al. (2016)
## https://doi.org/10.1177/0013164415584576
##
## authors: Timo
## data prepared by Uta & Sabine
## date: 03-07-2019
##
###############################################################################################
###############################################################################################

# Note: Includes several functionings for analyzing DTF that are used in the main analyses syntax.




# predicted response probabilities
IRT.p <- function(theta, xsi) {
  score <- IRT.score(xsi)
  d <- exp(score * theta - xsi)
  d / rowSums(d, na.rm = TRUE)  
}

# expected item scores
IRT.s <- function(theta, xsi) {
  score <- IRT.score(xsi)
  rowSums(score * IRT.p(theta, xsi), na.rm = TRUE)
}

# expected test scores
IRT.t <- function(theta, xsi) {
  sum(IRT.s(theta, xsi))
}

# item scoring 
IRT.score <- function(xsi) {
  k <- ncol(xsi)
  score <- matrix(seq_len(k) - 1, nrow = nrow(xsi), ncol = k, byrow = TRUE)
  if (ncol(xsi) >= 3) score[!is.na(xsi[, 3]), ] <- score[!is.na(xsi[, 3]), ] * 0.5
  score[is.na(xsi)] <- NA
  score
}
  
# maximal test score
IRT.ts <- function(xsi) {
  score <- IRT.score(xsi)
  sum(apply(score, 1, max, na.rm = TRUE))
}

# list of imputed plausible parameter values
IRT.ipars <- function(xsi, se, R = 100) {
  xsi[, 1] <- NA
  f <- !is.na(xsi)
  g <- mvrnorm(R, mu = xsi[f], Sigma = diag(se[f]))
  if (is.vector(g)) g <- matrix(g, nrow = 1)
  xsi.imp <- list()
  for (i in seq_len(R)) {
    x <- xsi
    x[, 1] <- 0
    x[f] <- g[i, ]
    xsi.imp[[i]] <- x
  }
  xsi.imp
}

# transform vector of xsi parameters to matrix of thresholds
# IRT.xsi2mat <- function(xsi) {
#   items <- unique(substr(names(xsi), 0, 10))
#   k <- max(as.numeric(substr(names(xsi), 15, 15))) + 1
#   mat <- matrix(0, nrow = length(items), ncol = k)
#   for (i in seq_len(k - 1)) {
#     mat[, i + 1] <- xsi[paste0(items, "_Cat", i)]
#   }
#   t(apply(mat, 1, cumsum))
# }

# expected test scores for imputed item parameters
IRT.tmat <- function(xsi.list, se.list, theta, R = 100) {
  k <- length(xsi.list)
  ts <- list()
  for(i in seq_len(k)) {
    ipars <- IRT.ipars(xsi = xsi.list[[i]], se = se.list[[i]], R = R)
    ts[[i]] <- matrix(NA, nrow = length(theta), ncol = R)
    for (j in seq_len(R)) {
      ts[[i]][, j] <- sapply(theta, IRT.t, xsi = ipars[[j]])
    }
    rownames(ts[[i]]) <- theta
  }
  names(ts) <- paste0("gr", seq_len(k))
  ts
}

# estimate differential test functioning
IRT.dtf <- function(xsi.list, se.list, n = 21, R = 100) {
  
  # quadrature nodes and weights  
  gquad <- gauss.quad(n)
  q <- 0.5 * ((6 - -6) * gquad$nodes + -6 + 6)
  w <- -0.5 * (-6 - 6) * gquad$weights / 12 # scale to sum of 1
  
  # expected test scores
  TSmat <- IRT.tmat(xsi.list = xsi.list, se.list = se.list, theta = q, R = R)
  
  # signed and unsigned DTF for each theta
  sDTF <- uDTF <- list()
  for (i in seq_len(length(TSmat) - 1)) {
    for (j in seq(i + 1, length(TSmat))) {
      if (i == j) next
      sDTF[[paste0("grp", j, i)]] <- uDTF[[paste0("grp", j, i)]] <- matrix(NA, nrow = n, ncol = R)
      for (r in seq_len(R)) {
        sDTF[[paste0("grp", j, i)]][, r] <- TSmat[[paste0("gr", j)]][, r] - TSmat[[paste0("gr", i)]][, r]
        uDTF[[paste0("grp", j, i)]][, r] <- abs(TSmat[[paste0("gr", j)]][, r] - TSmat[[paste0("gr", i)]][, r])
      }
    }
  }
  
  # signed and unsigend DTF across thetas
  wmat <- matrix(w, nrow = n, ncol = R, byrow = FALSE)
  DTF <- list()
  for (i in seq_len(length(sDTF))) {
    x <-  rbind(colSums(wmat * sDTF[[i]]), colSums(wmat * uDTF[[i]]))
    rownames(x) <- c("sDTF", "uDTF")
    DTF[[i]] <- x
  }
  names(DTF) <- names(sDTF)
  
  # return as parameter and standard error
  out <- list()
  out$sDTF <- t(sapply(DTF, function(x) { c(b = mean(x["sDTF", ]), se = sd(x["sDTF", ]))  }))
  out$uDTF <- t(sapply(DTF, function(x) { c(b = mean(x["uDTF", ]), se = sd(x["uDTF", ]))  }))
  out$sDTF21 <- cbind(b = rowMeans(sDTF$grp21), se = apply(sDTF$grp21, 1, sd))
  out$sDTF31 <- cbind(b = rowMeans(sDTF$grp31), se = apply(sDTF$grp31, 1, sd))
  out$sDTF41 <- cbind(b = rowMeans(sDTF$grp41), se = apply(sDTF$grp41, 1, sd))
  out$sDTF32 <- cbind(b = rowMeans(sDTF$grp32), se = apply(sDTF$grp32, 1, sd))
  out$sDTF42 <- cbind(b = rowMeans(sDTF$grp42), se = apply(sDTF$grp43, 1, sd))
  out$sDTF43 <- cbind(b = rowMeans(sDTF$grp43), se = apply(sDTF$grp43, 1, sd))
  out$uDTF21 <- cbind(b = rowMeans(uDTF$grp21), se = apply(uDTF$grp21, 1, sd))
  out$uDTF31 <- cbind(b = rowMeans(uDTF$grp31), se = apply(uDTF$grp31, 1, sd))
  out$uDTF41 <- cbind(b = rowMeans(uDTF$grp41), se = apply(uDTF$grp41, 1, sd))
  out$uDTF32 <- cbind(b = rowMeans(uDTF$grp32), se = apply(uDTF$grp32, 1, sd))
  out$uDTF42 <- cbind(b = rowMeans(uDTF$grp42), se = apply(uDTF$grp43, 1, sd))
  out$uDTF43 <- cbind(b = rowMeans(uDTF$grp43), se = apply(uDTF$grp43, 1, sd))
  rownames(out$sDTF21) <- rownames(out$sDTF31) <- rownames(out$sDTF41) <- rownames(out$sDTF32) <- 
    rownames(out$sDTF42) <- rownames(out$sDTF43) <- rownames(out$uDTF21) <- rownames(out$uDTF31) <- 
    rownames(out$uDTF41) <- rownames(out$uDTF32) <- rownames(out$uDTF42) <- rownames(out$uDTF43) <- 
    round(q, 1)
  out
}


# multi-facetted Rasch analyses for two groups
# returns DIF and main effects
IRT.dif <- function(resp, group, wgt) {
  
  # facets
  facets <- data.frame(group = unclass(group))

  if (max(resp, na.rm = TRUE) > 1) {
    formulaA = ~ item + item:step + group + item:group
    des <- designMatrices.mfr2(resp = resp, facets = facets, formulaA = formulaA)
    resp2 <- des$gresp$gresp.noStep
    A <- des$A$A.3d[ , , -des$xsi.elim[, 2]]
    # 0.5 scoring
    B <- des$B$B.3d
    s <- (substring(rownames(B), 8, 8) == "s") # PCMs
    B[s, , 1] <- B[s, , 1] / 2
  } else {
    formulaA = ~ item + group + item:group
    des <- designMatrices.mfr(resp = resp, facets = facets, formulaA = formulaA)
    resp2 <- des$gresp$gresp.noStep
    A <- des$A$A.3d
    B <- des$B$B.3d
  }
  
  # fit model
  fit <- tam.mml(resp = resp2, A = A, B = B, control = list(increment.factor = 1.02), 
                 pweight = wgt, verbose = FALSE)
  
  # differences in item parameters (positive values indiciate more difficult items in group 1 as compared to group 2)
  xsi <- data.frame(b12 = fit$xsi[grepl(":group1", rownames(fit$xsi)), 1] - # difference between groups 1 and 2
                          fit$xsi[grepl(":group2", rownames(fit$xsi)), 1])
  xsi$se12 <- (fit$xsi[grepl(":group2", rownames(fit$xsi)), 2]^2 +          # pooled standard error for groups 1 and 2
               fit$xsi[grepl(":group1", rownames(fit$xsi)), 2]^2)^0.5
  xsi$b13 <- fit$xsi[grepl(":group1", rownames(fit$xsi)), 1] -              # difference between groups 1 and 3
             fit$xsi[grepl(":group3", rownames(fit$xsi)), 1]
  xsi$se13 <- (fit$xsi[grepl(":group3", rownames(fit$xsi)), 2]^2 +          # pooled standard error for groups 1 and 3
               fit$xsi[grepl(":group1", rownames(fit$xsi)), 2]^2)^0.5
  xsi$b14 <- 2 * fit$xsi[grepl(":group1", rownames(fit$xsi)), 1]            # difference between groups 1 and 4
  xsi$se14 <- (fit$xsi[grepl(":group1", rownames(fit$xsi)), 2]^2 * 2)^0.5   # pooled standard error for groups 1 and 4
  xsi$b23 <- fit$xsi[grepl(":group2", rownames(fit$xsi)), 1] -              # difference between groups 2 and 3
             fit$xsi[grepl(":group3", rownames(fit$xsi)), 1]
  xsi$se23 <- (fit$xsi[grepl(":group3", rownames(fit$xsi)), 2]^2 +          # pooled standard error for groups 2 and 3
               fit$xsi[grepl(":group2", rownames(fit$xsi)), 2]^2)^0.5
  xsi$b24 <- 2 * fit$xsi[grepl(":group2", rownames(fit$xsi)), 1]            # difference between groups 2 and 4
  xsi$se24 <- (fit$xsi[grepl(":group2", rownames(fit$xsi)), 2]^2 * 2)^0.5   # pooled standard error for groups 2 and 4
  xsi$b34 <- 2 * fit$xsi[grepl(":group3", rownames(fit$xsi)), 1]            # difference between groups 3 and 4
  xsi$se34 <- (fit$xsi[grepl(":group3", rownames(fit$xsi)), 2]^2 * 2)^0.5   # pooled standard error for groups 3 and 4
  f <- rownames(fit$xsi)[grepl("_c$", rownames(fit$xsi))]
  rownames(xsi) <- f[-length(f)]
  
  # main effects (positive values indicate higher mean in group1 as compared to group 2)
  me <- data.frame(b12 = (fit$xsi["group2", 1] - fit$xsi["group1", 1])) # main effects for groups 1 and 2
  me$se12 <- (fit$xsi["group2", 2]^2 + fit$xsi["group1", 2]^2)^0.5      # pooled standard error for groups 1 and 2
  me$b13 <- (fit$xsi["group3", 1] - fit$xsi["group1", 1])               # main effect for groups 1 and 3
  me$se13 <- (fit$xsi["group3", 2]^2 + fit$xsi["group1", 2]^2)^0.5      # pooled standard error for groups 1 and 3
  me$b14 <- (2 * fit$xsi["group1", 1]) * -1                             # main effect for groups 1 and 4
  me$se14 <- (2 * fit$xsi["group1", 2]^2)^0.5                           # pooled standard error for groups 1 and 4
  me$b23 <- (fit$xsi["group3", 1] - fit$xsi["group2", 1])               # main effect for groups 2 and 3
  me$se23 <- (fit$xsi["group3", 2]^2 + fit$xsi["group2", 2]^2)^0.5      # pooled standard error for groups 2 and 3
  me$b24 <- (2 * fit$xsi["group2", 1]) * -1                             # main effect for groups 2 and 4
  me$se24 <- (2 * fit$xsi["group2", 2]^2)^0.5                           # pooled standard error for groups 2 and 4
  me$b34 <- (2 * fit$xsi["group3", 1]) * -1                             # main effect for groups 1 and 3
  me$se34 <- (2 * fit$xsi["group3", 2]^2)^0.5                           # pooled standard error for groups 1 and 2
  rownames(me) <- "main effect"
  
  rbind(xsi, me)
}
