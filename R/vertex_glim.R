mni.read.glim.file <- function(filename, header=FALSE, fill=FALSE) {
  glim <- as.data.frame(read.table(filename, header=header, fill=fill))
  return(glim)
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
  
# run stats at every vertex
mni.vertex.statistics <- function(glim.matrix, statistics.model=NA,
                                  vertex.table=FALSE) {
  # number of rows in the matrix
  number.subjects <- nrow(glim.matrix)
  # number of vertices. Assume that they are all the same, so just read it
  # off the first subject.
  number.vertices <-
    length(as.matrix(read.table(as.character(glim.matrix[1,1]))))


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
                    ncol = number.terms)),
                  q.values = data.frame(matrix(NA, nrow = number.vertices,
                    ncol = number.terms)))

  # assign the correct names
  names(results$slope) <- variable.names
  names(results$std.error) <- variable.names
  names(results$tstatistic) <- variable.names
  names(results$q.values) <- variable.names

  modulo <- 1000
  f <- formula(statistics.model)

  # run the stats at each vertex
  cat("   Percent done: ")
  for (i in 1:number.vertices) {
    y <- vertex.table[i,]
    s = summary(lm(f))
    results$adj.r.squared[i] <- s$adj.r.squared
    results$fstatistic[i] <- s$fstatistic[1]
    results$intercept[i] <- s$coefficients[1,1]
    for (j in 1:number.terms) {
      results$slope[i,j] <- s$coefficients[j,1]
      results$std.error[i,j] <- s$coefficients[j,2]
      results$tstatistic[i,j] <- s$coefficients[j,3]
    }
    # print progress report to the terminal
    if (i %% modulo == 0) {
      cat(format((i/number.vertices)*100, digits=3))
      cat("%  ")
    }
  }
  cat("\n")

  # compute the q values for all of the corresponding t-stats
  cat("   Computing q values\n")
  for (i in 1:number.terms) {
    q <- mni.compute.FDR(t.stats=results$tstatistic[,i],
                         df=number.subjects-1)
    results$q.values[,i] <- q$q
  }
  
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
  
  
