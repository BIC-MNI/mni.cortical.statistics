read.glim.file <- function(filename, header=FALSE) {
  glim <- as.data.frame(read.table(filename, header=header))
  return(glim)
}

# run statistics on the mean of all the files included in the glim model.
mean.statistics <- function(glim.matrix, statistics.model=NA) {
  # number of rows in the matrix
  l <- length(glim.matrix[,1])
  # build the table to hold all the means
  means <- matrix(NA, nrow = l, ncol = 1)
  for (i in 1:l) {
    means[i] <- mean(as.matrix(read.table(as.character(glim.matrix[i,1]))))
  }
  attach(glim.matrix)
  y <- means
  result <- lm(formula(statistics.model))
  return(result)
}
  
# run stats at every vertex
vertex.statistics <- function(glim.matrix, statistics.model=NA) {
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

  # get the number of terms in the formula
  # NOTE: assumes that terms are separated by a +, which won't always be true
  number.terms <- length(strsplit(statistics.model, '\\+')[[1]])
  
  # create the output table
  results <- list(adj.r.squared = vector(length = number.vertices),
                  fstatistic = vector(length = number.vertices),
                  intercept = vector(length = number.vertices),
                  slope = matrix(NA, nrow = number.vertices,
                    ncol = number.terms-1),
                  tstatistic = matrix(NA, nrow = number.vertices,
                    ncol = number.terms-1))

  # attach the named variables
  attach(glim.matrix)
  
  # run the stats at each vertex
  for (i in 1:number.vertices) {
    y <- vertex.table[i,]
    s = summary(lm(formula(statistics.model)))
    results$adj.r.squared[i] <- s$adj.r.squared
    results$fstatistic[i] <- s$fstatistic[1]
    results$intercept[i] <- s$coefficients[1,1]
    for (j in 2:number.terms) {
      results$slope[i,j-1] <- s$coefficients[j,1]
      results$tstatistic[i,j-1] <- s$coefficients[j,3]
    }
    print(i)
  }
  
  return(results)
  
}

# write out the statistics to a flat text file
write.vertex.stats <- function(vertex.stats, filename, headers = TRUE,
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
  


  
  
