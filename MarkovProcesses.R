##########################################################################
###########     Markov Process Simulation    #############################
##########################################################################
## General setup is a generating matrix A on a number of states n
## A number of possible models can be created including:
##        Poisson processes     - generatePPMatrix(n, lambda)
##        Queues on K servers   - generateQMatrix(n, K, lambda, mu)
##        Birth-death processes - generateBDMatrix(n, birthrate, deathrate)
###########################################################################
## Functions for simulation:
##      simulateMProcess(A, T, initial.state = 1, rounding = 6) 
##                        where A is the generating matrix,
##                        T is the maximum duration of the simulation
##                        output is a data frame with time/state variables
##                        
##      addnoise(x, error.time, error.state) adds noise to the data for 
##                                           authenticity (don't use!)
##
##      etimes(x, n) produces a data frame of empirical transition times
##      etransitionmatrix(x, n) produces the empirical transition times
##########################################################################
##########################################################################


### generates a birth-death matrix. The functions must return the birth/death rate for index n
generateBDMatrix <- function(n, birthrate, deathrate){
  if (!is.function(birthrate) | !is.function(deathrate)) {stop("birthrate and deathrate should be functions")}
  M <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    if (i > 1) { M[i,i-1] <- deathrate(i-1) }
    if (i < n) { M[i,i+1] <- birthrate(i-1) }
    M[i,i] <- -sum(M[i,])
  }
  return(M)
}

transitionMatrix <- function(A) {
  if (!is.matrix(A)) {stop("Input in transitionMatrix must be a square matrix")}
  n <- sqrt(length(A))
  if (n != round(n)) {stop("Input in transitionMatrix must be a square matrix")}
  P <- matrix(seq(length.out = n^2),nrow = n)
  for (i in 1:n) {
    denominator <- -A[i,i]
    for (j in 1:n) {
      if (j == i) { P[i,j] <- 0 } else { P[i,j] <- A[i,j]/denominator}
    }
  }
  return(P)
}

simulateMProcess <- function(A, T, initial.state = 1, rounding = 6) {
  ### A is a generator matrix, s is the initial state, T is the duration
  ### The function outputs a data frame with two variables "times" and "states"
  ### The rounding returns a rounded value for the times
  P <- transitionMatrix(A)
  n <- sqrt(length(A))
  lambda <- c()
  for (i in 1:n) {lambda[i] <- -A[i,i]}
  t <- 0
  state <- initial.state
  times <- c(t)
  states <- c(initial.state)
  while (t <= T) {
    t <- t + rexp(1,rate = lambda[state])
    times = append(times, round(t,rounding))
    
    transprob <- runif(1)
    for (k in 1:n) {
      if ( transprob < P[state,k]) {
        state <- k
        break
      }
      transprob <- transprob - P[state,k]
    }
    states <- append(states, state)
  }
  return(data.frame(times, states))
}

etransitionMatrix <- function(x,n) { #calculates the empirical transition matrix
  T <- matrix(0, nrow = n, ncol = n)
  previousstate <- x$states[1]
  for (i in 2:(length(x$times))) {
    currentstate <- x$states[i]
    T[previousstate,currentstate] <- T[previousstate,currentstate] + 1
    previousstate <- currentstate
  }
  P <- matrix(0,nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      P[i,j] <- T[i,j]/sum(T[i,])
    }
  }
return(P);
}

etimes <- function(x, n) { # calculates the empirical times
  t <- matrix(nrow = n,ncol=(length(x$times)))
  len <- c(length.out = n)
  for (i in 1:n) { len[i] <- 1 }
  for (i in 1:(length(x$times)-1)) {
  t[x$states[i],len[x$states[i]]] <- x$times[i+1] - x$times[i]
  len[x$states[i]] <- len[x$states[i]] + 1
  }
  return(t)
}

generateQMatrix <- function(n, K, lambda, mu) { 
  # generates a matrix for a K-server queue up to size nxn
  M <- matrix(0,nrow = n, ncol = n)
  M[1,1] <- -lambda
  M[1,2] <- lambda
  if (K >= 2 & K <= n) {
    for (i in 2:K) {
      M[i,i-1] <- (i-1)*mu
      M[i,i] <- -(lambda + (i-1)*mu)
      M[i,i+1] <- lambda
    }
  }
  if (K <= n) {
    for (i in (K+1):(n-1)) {
      M[i,i-1] <- K*mu
      M[i,i] <- -(lambda + K*mu)
      M[i,i+1] <- lambda      
    }
    M[n,n-1] <- K*mu
    M[n,n] <- -K*mu
  }
  return(M)
}

generatePPMatrix <- function(n, lambda) {
  # generates a nxn matrix for a Poisson process with intensity lambda
  br <- function(n) {return(lambda)}
  dr <- function(n) {return(0)}
  return(generateBDMatrix(n, br, dr))
}

addnoise <- function(x, error.time, error.state) {
  # this function adds noise to the Markov process simulation (don't use it for student data)
  if (error.time > 0) {
    for (i in 2:length(x$times)) { x$times[i] <- max(x$times[i-1], x$times[i] + rnorm(1,0,error.time/3)) }
  }
  if (error.state > 0) { # states must be integer so just use + or - for now
    warning("Adding noise to states may violate assumptions of a simulated Markov process")
    for (i in 1:length(x$states)) { x$states[i] <- x$states[i] + sample(-1:1,1) }
  }
  return(x)
}
