mni.read.glim.file <- function(filename, header=FALSE, fill=FALSE,
                               file.type="space") {
  glim <- NULL
  if(file.type == "space") {
    glim <- as.data.frame(read.table(filename, header=header, fill=fill))
  }
  else if (file.type == "csv") {
    glim <- as.data.frame(read.csv(filename, header=header, fill=fill))
  }
  else {
    stop("File type must be either space or csv")
  }
  return(glim)
}

# run a normal model from beginning to end
mni.statistical.run <- function(input, model, output, header=FALSE, fill=FALSE) {
  gf <- mni.read.glim.file(input, header, fill)
  dt <- mni.build.data.table(gf)
  ms <- mni.mean.statistics(gf, model, dt)
  vs <- mni.vertex.statistics(gf, model, dt)
  mni.write.vertex.stats(vs, output, mean.stats=ms, glim.matrix=gf)
  #return(list(gf=gf, ms=ms, vs=vs, dt=dt))
}

# run statistics on the mean of all the files included in the glim model.
mni.mean.statistics <- function(glim.matrix, statistics.model=NA,
                                vertex.table=FALSE) {
  # number of rows in the matrix
  l <- length(glim.matrix[,1])
  # build the table to hold all the means
  means <- matrix(NA, nrow = l, ncol = 1)

  # build the table to hold all of the values - unless they are given as
  # an argument
  if (mode(vertex.table) == "logical") {
    vertex.table <- mni.build.data.table(glim.matrix)
  }
  for (i in 1:l) {
    means[i] <- mean(vertex.table[,i])
  }
  attach(glim.matrix)
  y <- means
  result <- lm(formula(statistics.model))
  return(result)
}

# do a power analysis at every vertex
mni.vertex.power.analysis <- function(std.file, n=25, alpha=0.005,
                                      power=0.995, delta=0.5) {
  # read the info - the standard deviation file.
  vertex.table <- as.matrix(read.table(as.character(std.file)))
  number.vertices <- length(vertex.table)

  # create a list of two vectors to hold the results
  results <- list(delta.at.n = vector(length = number.vertices),
                  n.at.delta = vector(length = number.vertices))

  modula <- 0
  for (i in 1:number.vertices) {
      # the power analysis fails if the std is too low.
      if(vertex.table[i] < 0.02) {
        results$n.at.delta[i] <- 0
        results$delta.at.n[i] <- 0
      }
      else {
        # find the necessary n at the set delta
        p <- power.t.test(n=NULL, delta=delta, sd=vertex.table[i],
                          sig.level=alpha, power=power)
        results$n.at.delta[i] <- p$n
        # find the possible delta at a set n
        p <- power.t.test(n=n, delta=NULL, sd=vertex.table[i],
                          sig.level=alpha, power=power)
        results$delta.at.n[i] <- p$delta
      }
      tmp <- i %/% 1000
      if (tmp > modula) {
        print((i / number.vertices) * 100)
        modula <- tmp
      }
    }

  return(results)
}

# build a table for all the entries in each file
mni.build.data.table <- function(glim.matrix) {
  # number of rows in the matrix
  number.subjects <- length(glim.matrix[,1])
  # number of vertices. Assume that they are all the same, so just read it
  # off the first subject.
  number.vertices <-
    length(as.matrix(read.table(as.character(glim.matrix[1,1]))))

  # build the table to hold all of the values
  vertex.table <- matrix(NA, nrow = number.vertices, ncol = number.subjects)
  for (i in 1:number.subjects) {
    vertex.table[,i] <- as.matrix(read.table(as.character(glim.matrix[i,1])))
  }
  return(vertex.table)
}

# compare two different models at each vertex using an anova
mni.vertex.compare.models <- function(glim.matrix, model.one,
                                      model.two, vertex.table) {
  number.subjects <- nrow(glim.matrix)
  number.vertices <- nrow(vertex.table)

  # for debugging only
  #number.vertices <- 10

  attach(glim.matrix)
  
  results <- list(fstatistic = vector(length=number.vertices),
                  )

  for (i in 1:number.vertices) {
    y <- vertex.table[i,]
    lm.one <- lm(formula(model.one))
    lm.two <- lm(formula(model.two))
    a <- anova(lm.one, lm.two)
    results$fstatistic[i] <- a$F[2]
  }
  return(results)
}

# run a mixed effects model at every vertex
mni.vertex.mixed.model <- function(glim.matrix, fixed.effect, random.effect,
                                   vertex.table=FALSE) {
  # build the table to hold all of the values - unless they are given as
  # an argument
  if (mode(vertex.table) == "logical") {
    vertex.table <- mni.build.data.table(glim.matrix)
  }

  number.vertices <- nrow(vertex.table)

  attach(glim.matrix)
  # get the number of terms
  y <- vertex.table[1,]
  l <- lme(formula(fixed.effect), random=formula(random.effect))
  s <- summary(l)$tTable

  number.terms <- nrow(s)

  variable.names <- rownames(s)
  # remove the parentheses around the intercept term, as it is ugly when
  # written to file
  variable.names <- gsub('\[\(\)]', '', variable.names, perl=TRUE)

  # construct the output matrices
  value <- matrix(data=0, nrow=number.vertices, ncol=number.terms)
  std.error <- matrix(data=0, nrow=number.vertices, ncol=number.terms)
  t.value <- matrix(data=0, nrow=number.vertices, ncol=number.terms)

  modulo <- 500
  fe <- formula(fixed.effect)
  re <- formula(random.effect)
  # run the model at each vertex
  cat("    Percent done: ")
  for (v in 1:number.vertices) {
    y <- vertex.table[v,]
    s = try(summary(lme(fe, random=re))$tTable)

    # catch errors and blythely ignore them
    if (!inherits(s, "try-error")) {
      value[v,] <- s[,1]
      std.error[v,] <- s[,2]
      t.value[v,] <- s[,4]
    }

    # print progress report to the terminal
    if (v %% modulo == 0) {
      cat(format((v/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")

  # assign the correct names
  colnames(value) <- variable.names
  colnames(std.error) <- variable.names
  colnames(t.value) <- variable.names

  # create the output frame
  results <- list(value=value, std.error=std.error, t.value=t.value)
  return(results)
}


# test for homoscedasticity (I love that word)
mni.vertex.homoscedasticity <- function(glim.matrix, model, grouping,
                                        vertex.table) {
  number.vertices <- nrow(vertex.table)
  l.ratio <- vector(length=number.vertices)

  attach(glim.matrix)

  modulo <- 100
  
  for (i in 1:number.vertices) {
    y <- vertex.table[i,]
    gls1 <- gls(formula(model))
    gls2 <- gls(formula(model), weights = varIdent(form = formula(grouping)))
    l.ratio[i] <- 2 * abs(diff(c(logLik(gls1), logLik(gls2))))
    
    # print progress report to the terminal
    if (i %% modulo == 0) {
      cat(format((i/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")

  return(l.ratio)
}


# compare two different models at each vertex using an anova
mni.vertex.mixed.model.compare.models <- function(glim.matrix,
                                                  model.one,
                                                  model.two,
                                                  random.effect,
                                                  vertex.table) {
  number.subjects <- nrow(glim.matrix)
  number.vertices <- nrow(vertex.table)

  # for debugging only
  #number.vertices <- 10

  attach(glim.matrix)
  
  results <- list(l.ratio = vector(length=number.vertices),
                  p.value = vector(length=number.vertices))
                  
  
  for (i in 1:number.vertices) {
    y <- vertex.table[i,]

    model.one <- formula(model.one)
    model.two <- formula(model.two)
    r.effect <- formula(random.effect)

      
    l1 <- try(lm.one <- lme(model.one, random=r.effect, method="ML"))
    l2 <- try(lm.two <- lme(model.two, random=r.effect, method="ML"))
    if (!inherits(l1, "try-error") && !inherits(l2, "try-error")) {
      lratio <- 2 * abs(diff(c(lm.one$logLik, lm.two$logLik)))
      results$l.ratio[i] <- lratio
      results$p.value[i] <- 1 - pchisq(lratio, lm.one$fixDF$X[2] - lm.two$fixDF$X[2])
    }
    else {
      results$l.ratio[i] <- 0
      results$p.value[i] <- 1
    }
    
#    lm.one <- lme(y ~ Age, random=r.effect, method="ML")
#    lm.two <- lme(y ~ Age + I(Age^2), random=r.effect, method="ML")
#    a <- anova.lme(lm.one, lm.two)
#    results$l.ratio[i] <- a[2,8]
#    results$p.value[i] <- a[2,9]

    
  }
  return(results)
}


mni.vertex.mixed.model.anova <- function(glim.matrix, fixed.effect,
                                         random.effect, vertex.table=FALSE) {
  # build the table to hold all of the values - unless they are given as
  # an argument
  if (mode(vertex.table) == "logical") {
    vertex.table <- mni.build.data.table(glim.matrix)
  }

  number.vertices <- nrow(vertex.table)

  attach(glim.matrix)
  # get the number of terms
  y <- vertex.table[1,]
  l <- lme(formula(fixed.effect), random=formula(random.effect))
  a <- anova(l)

  number.terms <- nrow(a)

  variable.names <- rownames(a)
  # remove the parentheses around the intercept term, as it is ugly when
  # written to file
  variable.names <- gsub('\[\(\)]', '', variable.names, perl=TRUE)

  # construct the output matrices
  f.stats <- matrix(data=0, nrow=number.vertices, ncol=number.terms)

  modulo <- 500
  fe <- formula(fixed.effect)
  re <- formula(random.effect)
  # run the model at each vertex
  cat("    Percent done: ")
  for (v in 1:number.vertices) {
    y <- vertex.table[v,]
    a = try(anova(lme(fe, random=re))$"F-value")

    # catch errors and blythely ignore them
    if (!inherits(a, "try-error")) {
      f.stats[v,] <- a
    }

    # print progress report to the terminal
    if (v %% modulo == 0) {
      cat(format((v/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")

  # assign the correct names
  colnames(f.stats) <- variable.names

  # create the output frame
  results <- list(f.stats=f.stats)
  return(results)
}

# run an anova at every vertex
mni.vertex.anova <- function(glim.matrix, statistics.model=NA,
                             vertex.table=FALSE) {
  # build the table to hold all of the values - unless they are given as
  # an argument
  if (mode(vertex.table) == "logical") {
    vertex.table <- mni.build.data.table(glim.matrix)
  }

  number.subjects <- nrow(glim.matrix)
  number.vertices <- nrow(vertex.table)

  # attach the named variables
  attach(glim.matrix)

  # get the number of terms in the formula
  # run one anova.
  y <- vertex.table[1,]
  a <- aov(formula(statistics.model))
  s <- summary(a)

  variable.names <- rownames(s[[1]])
  # remove the residuals term
  variable.names <- variable.names[-(length(variable.names))]
  number.terms <- length(variable.names)
  
  #create the output table holding the F statistics
  results <- matrix(data=NA, nrow=number.vertices, ncol=number.terms)

  modulo <- 1000
  f <- formula(statistics.model)

  # run the anova at each vertex
  cat("    Percent done: ")
  for (v in 1:number.vertices) {
    y <- vertex.table[v,]
    s <- summary(aov(f))
    results[v,] <- s[[1]]$"F value"[1:number.terms]
    # print progress report to the terminal
    if (v %% modulo == 0) {
      cat(format((v/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")
  colnames(results) <- variable.names

  return(results)
}

# get the residuals of a linear model
mni.vertex.residuals <- function(glim.matrix, statistics.model, vertex.table) {

  number.vertices <- nrow(vertex.table)
  new.dt <- vertex.table
  attach(glim.matrix)
  modulo <- 1000
  # run the stats at each vertex
  cat("   Percent done: ")
  for (i in 1:number.vertices) {
    y <- vertex.table[i,]
    l <- lm(formula(statistics.model))
    new.dt[i,] <- residuals(l)
    # print progress report to the terminal
    if (i %% modulo == 0) {
      cat(format((i/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")

  detach(glim.matrix)

  return(new.dt)
}

# trace anatomical connectivity using several hops
mni.anatcon.trace <- function(data.table, y, cortex, hops=3, min.distance=40,
                              min.value=0.6) {

  filename.cor <- "/tmp/cor.vertstats"
  filename.peaks <- "/tmp/peaks.csv"
  
  cor1 <- mni.vertex.correlation(data.table, y)
  mni.write.vertex.stats(cor1, filename.cor);
  system(paste("vertstats_find_peaks -min_value", min.value,
               "-min_distance", min.distance, filename.cor,
               cortex, filename.peaks))
  peaks <- read.csv(filename.peaks, header=TRUE)

  cors <- matrix(ncol=length(peaks$vertex), nrow=length(cor1))
  for (i in 1:length(peaks$vertex)) {
    cat(paste("Correlating peak:", i, "\n"))
    cors[,i] <- mni.vertex.correlation(data.table, dt[peaks$vertex[i]+1,])
  }

  m <- apply(cors, 1, max)
  output <- ((cor1 > min.value) * 2) + ( m > min.value)
  
  return(output)
}

# get the strength of cross cortex correlations at every vertex
mni.vertex.correlation.strength <- function(data.table) {
  number.vertices <- nrow(data.table)

  modulo <- 10
  results <- vector(length=number.vertices)
  for (i in 1:number.vertices) {
    c <- mni.vertex.correlation(data.table, data.table[i,])
    results[i] <- sum(c)
    if (i %% modulo == 0) {
      cat(format((i/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")
  return(results)
}

# run a permutation test for significance of correlation
mni.vertex.correlation.permutation <- function(data.table, y, groups, num=100) {

  number.subjects <- ncol(data.table)
  results <- vector(length=num)
  for (i in 1:num) {
    g <- sample(groups)
    c1 <- mni.vertex.correlation(data.table[,g==levels(g)[1]],
                                 y[g==levels(g)[1]])
    c2 <- mni.vertex.correlation(data.table[,g==levels(g)[3]],
                                 y[g==levels(g)[3]])
    results[i] <- min(c1-c2)

    cat("Permutation: ")
    cat(i)
    cat("\n")
  }
  return(results)
}
  

# correlate every vertex with variable y
mni.vertex.correlation <- function(data.table, y) {

  number.subjects <- ncol(data.table)
  number.vertices <- nrow(data.table)

  results <- vector(length=number.vertices)

  for (i in 1:number.vertices) {
    results[i] <- cor(data.table[i,], y)
  }
  return(results)
}

# a variance test between two subsets. Optionally computes the
# variance in a robust fashion using bootstrap sampling of the
# residuals. To enable this set the bootstrap argument to be the ratio
# of the data's length to be sampled each time.
mni.vertex.var.test <- function(glim.matrix, statistics.model, subset1,
                               subset2, vertex.table, bootstrap=NULL) {

  number.subjects <- nrow(glim.matrix)
  number.vertices <- nrow(vertex.table)
  variance.ratio <- vector(length=number.vertices)
  f.stat <- vector(length=number.vertices)
  p.value <- vector(length=number.vertices)

  attach(glim.matrix)
  modulo <- 1000

  if(! is.null(bootstrap)) { # some variables for robust variance estimation
    m <- 100
    res <- numeric(m)
  }

  for (i in 1:number.vertices) {
    y <- vertex.table[i,]
    l1 <- lm(formula(statistics.model), glim.matrix, subset=subset1)
    l2 <- lm(formula(statistics.model), glim.matrix, subset=subset2)

    if (is.null(bootstrap)) { #default estimation of variance
      v <- var.test(l1, l2)
      variance.ratio[i] <- v$estimate
      f.stat[i] <- v$statistic
      p.value[i] <- v$p.value
    }
    else { # more robust estimation of variance
      res1 <- resid(l1)
      res2 <- resid(l2)
      for (v in 1:m) res[v] <- var(sample(res1,
                                          size=round(length(res1 / bootstrap)),
                                          replace=T))
      var1 <- median(res)
      for (v in 1:m) res[v] <- var(sample(res2,
                                          size=round(length(res2 / bootstrap)),
                                          replace=T))
      var2 <- median(res)
      vr <- var1 / var2
      variance.ratio[i] <- vr
      f.stat[i] <- vr
      p.value[i] <- 0
    }

    # print progress report to the terminal
    if (i %% modulo == 0) {
      cat(format((i/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")
  detach(glim.matrix)
  return(data.frame(f.stat, variance.ratio, p.value))
}

mni.vertex.mood.test <- function(glim.matrix, statistics.model, subset1,
                                 subset2, vertex.table) {

  number.subjects <- nrow(glim.matrix)
  number.vertices <- nrow(vertex.table)
  #variance.ratio <- vector(length=number.vertices)
  z.stat <- vector(length=number.vertices)
  p.value <- vector(length=number.vertices)

  attach(glim.matrix)
  modulo <- 1000

  for (i in 1:number.vertices) {
    y <- vertex.table[i,]
    l1 <- lm(formula(statistics.model), glim.matrix, subset=subset1)
    l2 <- lm(formula(statistics.model), glim.matrix, subset=subset2)
    v <- mood.test(residuals(l1), residuals(l2))
    #variance.ratio[i] <- v$estimate
    z.stat[i] <- v$statistic
    p.value[i] <- v$p.value
    # print progress report to the terminal
    if (i %% modulo == 0) {
      cat(format((i/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")
  detach(glim.matrix)
  return(data.frame(z.stat,p.value))
}


# run stats at every vertex
mni.vertex.statistics <- function(glim.matrix, statistics.model=NA,
                                  vertex.table=FALSE) {

  # number of rows in the matrix
  number.subjects <- nrow(glim.matrix)
  # number of vertices. Assume that they are all the same, so just read it
  # off the first subject.
  if (mode(vertex.table) != "logical") {
    number.vertices <- nrow(vertex.table)
  }
  else {
    number.vertices <-
      length(as.matrix(read.table(as.character(glim.matrix[1,1]))))
  }

  # build the table to hold all of the values - unless they are given as
  # an argument
  if (mode(vertex.table) == "logical") {
    vertex.table <- matrix(NA, nrow = number.vertices, ncol = number.subjects)
    for (i in 1:number.subjects) {
      vertex.table[,i] <- as.matrix(read.table(as.character(glim.matrix[i,1])))
    }
  }

  #number.vertices <- 500
  
  # attach the named variables
  attach(glim.matrix)

  # get the number of terms in the formula
  # run one lm to determine all the required info
  y <- vertex.table[1,]
  l <- lm(formula(statistics.model))
  s <- summary(l)
  variable.names <- row.names(s$coefficients)
  # remove the parentheses around the intercept term, as it is ugly when
  # written to file
  variable.names <- gsub('\[\(\)]', '', variable.names, perl=TRUE)
  number.terms <- length(variable.names)

  # create the output table
  results <- list(adj.r.squared = vector(length = number.vertices),
                  fstatistic = vector(length = number.vertices),
                  intercept = vector(length = number.vertices),
                  slope = data.frame(matrix(NA, nrow = number.vertices,
                    ncol = number.terms)),
                  std.error = data.frame(matrix(NA, nrow = number.vertices,
                    ncol = number.terms)),
                  tstatistic = data.frame(matrix(NA, nrow = number.vertices,
                    ncol = number.terms)))
   

  slope <- matrix(data=NA, nrow=number.vertices, ncol=number.terms)
  tstats <- matrix(data=NA, nrow=number.vertices, ncol=number.terms)
  stderr <- matrix(data=NA, nrow=number.vertices, ncol=number.terms)

  modulo <- 1000
  f <- formula(statistics.model)

  # run the stats at each vertex
  cat("   Percent done: ")
  for (v in 1:number.vertices) {
    y <- vertex.table[v,]
    s = summary(lm(f))
    results$adj.r.squared[v] <- s$adj.r.squared
    results$fstatistic[v] <- s$fstatistic[1]
    results$intercept[v] <- s$coefficients[1,1]

    slope[v,] <- s$coefficients[,1]
    stderr[v,] <- s$coefficients[,2]
    tstats[v,] <- s$coefficients[,3]
    
    # print progress report to the terminal
    if (v %% modulo == 0) {
      cat(format((v/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")

  results$slope = slope
  results$std.error = stderr
  results$tstatistic = tstats

  # assign the correct names
  colnames(results$slope) <- variable.names
  colnames(results$std.error) <- variable.names
  colnames(results$tstatistic) <- variable.names

  # compute the q values for all of the corresponding t-stats
  #cat("   Computing q values\n")
  #for (i in 1:number.terms) {
  #  q <- mni.compute.FDR(t.stats=results$tstatistic[,i],
  #                       df=number.subjects-1)
  #  results$q.values[,i] <- q$q
  #}

  return(results)

}

# write out the statistics to a flat text file
mni.write.vertex.stats <- function(vertex.stats, filename, headers = TRUE,
                               mean.stats = NULL, glim.matrix = NULL) {

  # write the header
  append.file = TRUE
  if (headers == TRUE) {
    # always clobber the file
    write("<header>", file = filename)
    
    if (is.object(mean.stats)) {
      write("<mean>", file = filename, append = TRUE)
      sink(filename, append = TRUE)
      print(summary(mean.stats))
      sink(NULL)
      write("</mean>", file = filename, append = TRUE)
      write("<formula>", file = filename, append = TRUE)
      sink(filename, append = TRUE)
      print(formula(mean.stats))
      sink(NULL)
      write("</formula>", file = filename, append = TRUE)
    }
    if (is.data.frame(glim.matrix)) {
      write("<matrix>", file = filename, append = TRUE)
      write.table(glim.matrix, file = filename, append = TRUE,
                  row.names = FALSE)
      write("</matrix>", file = filename, append = TRUE)
    }
    write("</header>", file = filename, append = TRUE)
  }
  else {
    append.file = FALSE
  }
  
  write.table(vertex.stats, file = filename, append = append.file,
              quote = FALSE, row.names = FALSE, col.names = headers)
}

# read a single vertstats column from file
mni.read.vertstats.column <- function(filename, column.name=NULL) {

  # by name
  if (is.null(column.name)) {
    # default name is equal to the first column
    column.name = "Column0"
  }

  # there is a bug in vertstats_extract, causing it to have one too many
  # values at the end - so this is a temporary fix. Ick!
  return.val <- as.numeric( system( paste("vertstats_extract", filename,
                                    column.name), intern=TRUE))
  l <- length(return.val)
  if (l == 40963) {
    return.val <- return.val[1:l-1]
  }
  return(return.val)
}

    
  
mni.compute.FDR <- function(t.stats=NULL, p.values=NULL, filename=NULL,
                            column.name=NULL,
                            df=Inf, fdr=0.05, plot.fdr=FALSE) {
  # argument handling: must have either t.stats or p.values
  if (is.null(t.stats) && is.null(p.values) && is.null(filename)) {
    stop("Either t.stats or p.values have to be specified")
  }
  if (! is.null(filename)) {
    if (is.null(column.name)) {
      data.headers <- gsub(' +', '', system(paste("vertstatsinfo -dataheaders",
                                                  filename), intern=TRUE),
                           perl=TRUE)
      cat(" Column Names of t-statistics: \n\n")
      data.headers <- data.headers[grep('tstatistic', data.headers)]
      print(data.headers)
      cat("\n")
      stop("Specify a column name with the file name, see choices above")
    }
    # get the actual stats
    t.stats <- as.numeric(system(paste("vertstats_extract", filename,
                                       column.name),
                                 intern=TRUE))
  }
                        
      
      
        
  if (is.null(p.values)) {
    # compute the p-values from the t-stats
    p.values <- abs(pt(abs(t.stats), df) - 1)
  }
  # sort the p-values
  sorted.p.values <- sort(p.values, index.return = TRUE)
  # compute the q stats
  q <- sorted.p.values$x / 1:length(p.values) * length(p.values)
  if (plot.fdr == TRUE) {
    plot(1:length(p.values)/length(p.values), sorted.p.values$x)
    abline(0, fdr, col="red")
  }
  # find the threshold
  q2 <- q <= fdr
  r <- sort(q2, decreasing = TRUE, index.return = TRUE)
  fdr.threshold <- qt(sorted.p.values$x[max(r$ix[r$x == TRUE])], df)
  # sort the q values to be in the same order as the t.stats passed in
  q[sorted.p.values$ix] <- q
  # return threshold and q values.
  return(fdr.threshold, q)
}
  
  
