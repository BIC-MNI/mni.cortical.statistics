
library(MASS)
library(nnet)

# check the ability of a single predictor to reproduce the true
# diagnosis using the leave one out method (i.e. it will do as many
# lda fits as there are subjects).
mni.leave.one.out.diagnosis <- function(values, true.diagnosis, method="lda") {
  number.cases <- length(true.diagnosis)
  diagnosis <- matrix(NA, ncol=1, nrow=number.cases)
  if (method == "qda") {
    diagnosis = qda(as.matrix(values), true.diagnosis, CV=TRUE)$class
  }
  else if (method == "lda") {
    diagnosis = lda(as.matrix(values), true.diagnosis, CV=TRUE)$class
  }
  else if (method == "multinom") {
    for (i in 1:number.cases) {
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
  

predplot.jpl <- function(x,y,truth, len=100, xlab="x", ylab="y") {
  r1 <- range(x)
  r2 <- range(y)
  s1 <- seq(r1[1], r1[2], length=len)
  s2 <- seq(r2[1], r2[2], length=len)

  grid <- expand.grid(x=s1, y=s2)

  l <- lda(cbind(x=x, y=y), truth)
  q <- qda(cbind(x=x, y=y), truth)
  m <- multinom(truth ~ x + y)
  
  pred.l <- predict(l, newdata=grid)
  pred.q <- predict(q, newdata=grid, type="predictive")
  pred.m <- predict(m, grid, type = "probs")

  plot(x,y,type="n", xlab=xlab, ylab=ylab)
  #text(x,y,labels=as.character(truth))
  points(x,y,col=as.numeric(truth), pch=19)
  contour(s1, s2, matrix(pred.l$post[,1], len), add=T,
          levels=0.5, drawlabels=F)
  contour(s1, s2, matrix(pred.q$post[,1], len), add=T, levels=0.5,
          drawlabels=F, col="red")
  contour(s1, s2, matrix(pred.m, len), add=T, levels=0.5,
          drawlabels=F, col="blue")

  legend(locator(n=1), c(levels(truth)[1], levels(truth)[2], "lda", "qda",
                    "logistic"),
         pch=c(19,19,-1, -1, -1), col=c(1, 2, "black", "red", "blue"),
         lty=c(0,0,1, 1,1))
  
}

predplot.single.jpl <- function(x, truth, len=100, xlab="x", ylab="y") {
  r1 <- range(x)
  s1 <- seq(r1[1], r1[2], length=len)

  n <- x
  
  l <- lda(as.matrix(x), truth)
  q <- qda(as.matrix(x), truth)
  m <- multinom(truth ~ n)

  n <- s1
  pred.l <- predict(l, newdata=as.matrix(s1))
  pred.q <- predict(q, newdata=as.matrix(s1), type="predictive")

  pred.m <- predict(m, newdata=as.matrix(n), type="probs")
  
  index.l <- which.min(abs(pred.l$post[,1] -0.5))
  index.q <- which.min(abs(pred.q$post[,1] -0.5))
  index.m <- which.min(abs(pred.m - 0.5))
  #plot(as.factor(truth), x)

  stripplot(x ~ as.factor(truth), xlab=xlab, ylab=ylab,
            panel=function(x,y) {
              panel.stripplot(x,y);
              panel.abline(h=s1[index.l]);
              panel.abline(h=s1[index.q], col="red");
              panel.abline(h=s1[index.m], col="blue")
            })

  #return(list(lda=s1[index.l], qda=s1[index.q], multinom=s1[index.m]))
}

simulate.discriminant <- function(group.means, group.sds, group.ns,
                                  num.sims=100, method="lda") {
  truth <- c(rep(1, group.ns[1]), rep(0, group.ns[2]))

  results <- list(sensitivity = vector(length=num.sims),
                  specificity = vector(length=num.sims),
                  accuracy = vector(length=num.sims))

  for (i in 1:num.sims) {
    data <- c(rnorm(group.ns[1], group.means[1], group.sds[1]),
              rnorm(group.ns[2], group.means[2], group.sds[2]))
    l <- lda(as.matrix(data), truth)
    t <- as.data.frame(table(predict(l)$class, truth))
    results$sensitivity[i] <- t[4,3] / (t[4,3] + t[3,3])
    results$specificity[i] <- t[1,3] / (t[1,3] + t[2,3])
    results$accuracy[i] <- (t[4,3] + t[1,3]) / (t[1,3] + t[2,3] + t[3,3] + t[4
,3])
  }
  return(results)
}

discriminate.on.areas <- function(data.table, animal.table, truth,
                                  method="lda") {
  areas <- unique(animal.table)
  num.areas <- length(areas)

  results <- list(area = vector(length=num.areas),
                  accuracy = vector(length=num.areas),
                  sensitivity = vector(length=num.areas),
                  specificity = vector(length=num.areas))
  
  for (i in 1:num.areas) {
    results$area[i] <- areas[i]
    if (sum(animal.table == areas[i]) > 1) {
      d <- mni.leave.one.out.diagnosis(colMeans(data.table[animal.table == areas[i],]),
                                       truth, method=method)
      r <- mni.vertex.sensitivity(t(d), truth)
      
      results$accuracy[i] <- r$accuracy
      results$sensitivity[i] <- r$sensitivity
      results$specificity[i] <- r$specificity
    }
  }
  return(results)
}

discriminate.two.areas <- function(data.table, animal.table, truth,
                                   method="lda") {
  areas <- unique(animal.table)
  combo <- expand.grid(areas, areas)
  num.combo <- nrow(combo)

  results <- list(area1 = vector(length=num.combo),
                  area2 = vector(length=num.combo),
                  accuracy = vector(length=num.combo),
                  sensitivity = vector(length=num.combo),
                  specificity = vector(length=num.combo),
                  accuracy.a1 = vector(length=num.combo),
                  accuracy.a2 = vector(length=num.combo))
  
  for (i in 1:num.combo) {
    results$area1[i] <- combo[i,1]
    results$area2[i] <- combo[i,2]
    print(i)
    if (combo[i,1] != combo[i,2] && sum(animal.table == combo[i,1]) > 1
        && sum(animal.table == combo[i,2]) > 1) {
      a1 <- colMeans(data.table[animal.table == combo[i, 1],])
      a2 <- colMeans(data.table[animal.table == combo[i, 2],])
      d <- mni.leave.one.out.diagnosis(cbind(a1,a2), truth, method=method)
      r <- mni.vertex.sensitivity(t(as.numeric(d)), truth)

      d1 <- mni.leave.one.out.diagnosis(a1, truth, method=method)
      d2 <- mni.leave.one.out.diagnosis(a2, truth, method=method)
      r1 <- mni.vertex.sensitivity(t(as.numeric(d1)), truth)
      r2 <- mni.vertex.sensitivity(t(as.numeric(d2)), truth)
      
      results$accuracy[i] <- r$accuracy
      results$sensitivity[i] <- r$sensitivity
      results$specificity[i] <- r$specificity
      results$accuracy.a1[i] <- r1$accuracy
      results$accuracy.a2[i] <- r2$accuracy
    }
  }
  return(results)
    
}
  
  
