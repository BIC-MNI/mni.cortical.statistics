
library(MASS)
library(nnet)

# check the ability of a single predictor to reproduce the true
# diagnosis using the leave one out method (i.e. it will do as many
# lda fits as there are subjects).
mni.leave.one.out.diagnosis <- function(values, true.diagnosis, method="lda") {
  number.cases <- length(true.diagnosis)
  diagnosis <- matrix(NA, ncol=1, nrow=number.cases)
  for (i in 1:number.cases) {
    if (method == "lda") {
      lda.fit <- lda(as.matrix(values[-c(i)]), true.diagnosis[-c(i)])
      diagnosis[i] <- predict(lda.fit, values[i])$class
    }
    else if (method == "qda") {
      qda.fit <- qda(as.matrix(values[-c(i)]), true.diagnosis[-c(i)])
      diagnosis[i] <- predict(qda.fit, values[i])$class
    }
    else if (method == "multinom") {
      multinom.fit <- multinom(true.diagnosis ~ values, subset=-i)
      tmp.diagnosis <- predict(multinom.fit, values)
      diagnosis[i] <- tmp.diagnosis[i]
    }
  }
  return(diagnosis)
}

# runs the mni.leave.one.out.diagnosis at each vertex using only that
# vertex as a predictor. Returns a matrix with nrows = nvertices and
# each column corresponding to the diagnosis of each subject at that
# vertex.
mni.diagnostic.capabilities.of.vertices <- function(data.table, true.diagnosis, method="lda") {
  number.cases <- length(true.diagnosis)
  number.vertices <- nrow(data.table)
  results <- matrix(NA, ncol=number.cases, nrow=number.vertices)
  for (i in 1:number.vertices) {
    try(results[i,] <- mni.leave.one.out.diagnosis(data.table[i,], true.diagnosis, method=method))
    print(i)
  }
  return(results)
}

mni.vertex.sensitivity <- function(diagnostic.results, true.diagnosis) {
  number.vertices = nrow(diagnostic.results)
  results <- list(sensitivity = vector(length=number.vertices),
                  specificity = vector(length=number.vertices),
                  ppv = vector(length=number.vertices),
                  npv = vector(length=number.vertices),
                  accuracy = vector(length=number.vertices))
  for (i in 1:number.vertices) {
    t <- try(as.data.frame(table(diagnostic.results[i,], true.diagnosis)))
    if (!inherits(t, "try-error")) {
      results$sensitivity[i] <- t[4,3] / (t[4,3] + t[3,3])
      results$specificity[i] <- t[1,3] / (t[1,3] + t[2,3])
      results$ppv[i] <- t[4,3] / (t[4,3] + t[2,3])
      results$npv[i] <- t[1,3] / (t[1,3] + t[3,3])
      results$accuracy[i] <- (t[4,3] + t[1,3]) / (t[1,3] + t[2,3] + t[3,3] + t[4,3])
    }
    print(i)
  }
  return(results)
}

mni.stepwise.diagnostics <- function(diagnostic.results, true.diagnosis,
                                     data.table, max.steps = 1000) {
  attach(diagnostic.results)
  combined <- sensitivity + specificity #+ ppv + npv
  sorted <- sort(combined, index.return=TRUE, decreasing=TRUE)
  i <- 2
  perfect <- FALSE
  outcomes <- matrix(NA, nrow=max.steps, ncol=4)
  colnames(outcomes) <- c("TN", "FN", "FP", "TP")
  while (i < max.steps && perfect == FALSE) {
    #l <- lda(t(data.table[sorted$ix[1:i],]), true.diagnosis)
    #pl <- predict(l, t(data.table[sorted$ix[1:i],]))
    diagnosis <- mni.leave.one.out.diagnosis(t(data.table[sorted$ix[1:i],]),
                                               true.diagnosis)
    results <- as.data.frame(table(diagnosis, true.diagnosis))
    if (results[2,3] == 0 && results[3,3] == 0) {
      perfect <- TRUE
    }
    outcomes[i,] <- results[,3]
    i <- i + 1
  }
  detach(diagnostic.results)
  return(i, outcomes)
}
  
