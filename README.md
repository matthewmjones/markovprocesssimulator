# Markov Process Simulation

## Functions for simulating the processes
- ```simulateMProcess(A, T, initial.state = 1, rounding = 6)``` where ```A``` is the generating matrix, ```T``` is the maximum duration of the simulation output is a data frame with time/state variables
- ```addnoise(x, error.time, error.state)``` adds noise to the data for authenticity (this often leads to violations of basic assumptions and is not recommended)
- ```etimes(x, n)``` produces a data frame of empirical transition times
- ```etransitionmatrix(x, n)``` produces the empirical transition times

## In build generating matrices
The general setup consists of a generating matrix ```A``` on a number of states ```n```. A number of possible generating matrices can be created using the following two functions:

- Poisson processes     - ```generatePPMatrix(n, lambda)```
- Queues on ```K``` servers   - ```generateQMatrix(n, K, lambda, mu)```
  
  Each of these uses the following general function
- Birth-death processes - ```generateBDMatrix(n, birthrate, deathrate)```

Here ```birthrate``` and ```deathrate``` are functions of the form:

```function(i) { <return birth/death rate for row i> }```


