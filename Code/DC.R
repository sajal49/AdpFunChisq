# The original codes for DC was obtained from:
# http://eda.mmci.uni-saarland.de/prj/cisc/

require(Rfast)
# for testing
require(assertthat)

dc <- function(X, Y){
  
  # calculate frequencies
  freq_x = table(X)
  freq_y = table(Y)
  freq_xy = table(X, Y)
  
  # calculate marginal probabilities
  prob_x = freq_x/sum(freq_x)
  prob_y = freq_y/sum(freq_y)
  
  # calculate cell probabilities
  cell_xtoy = freq_xy
  for(i in 1:length(freq_x))
    cell_xtoy[i,] = cell_xtoy[i,]/freq_x[i]
  cell_xtoy = t(cell_xtoy)
  cell_ytox = freq_xy
  for(j in 1:length(freq_y))
    cell_ytox[,j] = cell_ytox[,j]/freq_y[j]
  
  # get dcor values
  dXtoY = Rfast::dcor(x = prob_x, y = cell_ytox)
  dYtoX = Rfast::dcor(x = prob_y, y = cell_xtoy)
  
  # return directional stats
  if(is.na(dXtoY$dcor))
    dXtoY$dcor=0
  if(is.na(dYtoX$dcor))
    dYtoX$dcor=0
  return(list(x.to.y = dXtoY$dcor, y.to.x = dYtoX$dcor))
  
}


# 3 test cases with truth obtained from the dc.py
test_dc <- function(){
  
  # Test 1
  X = c(23, 23, 24, 2, 2, 3, 3, 34, 23)
  Y = c(1, 1, 1, 1, 0, 0, 0, 0, 1)
  dc_res = dc(X,Y)
  assert_that(are_equal(dc_res$x.to.y, 0.7714911689255104))
  assert_that(are_equal(dc_res$y.to.x, 1))
  
  # Test 2
  X = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)
  Y = c(0, 0, 0, 1, 1, 1, 0, 0, 1, 0)
  dc_res = dc(X,Y)
  assert_that(are_equal(dc_res$x.to.y, 0.7888280382727337))
  assert_that(are_equal(dc_res$y.to.x, 1))
  
  # Test 3
  X = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  Y = c(0, 0, 0, 1, 1, 1, 2, 2, 2)
  dc_res = dc(X,Y)
  assert_that(are_equal(dc_res$x.to.y, 0))
  assert_that(are_equal(dc_res$y.to.x, 0))
  
  # Test 4
  X = c(1,1,1, 2,2,2, 3,3,3, 1,1,1, 4,4,4,4)
  Y = c(2,2,2, 1,1,1, 2,2,2, 2,2,2, 1,1,1,1)
  dc_res = dc(X, Y)
  assert_that(are_equal(dc_res$x.to.y, 0.7340934979782122))
  assert_that(are_equal(dc_res$y.to.x, 1))
  
  # Test 5
  X = c(1,2,3, 3,2,1, 1,1,1, 2,2,2, 3,3,3)
  Y = c(1,1,1, 2,2,2, 1,2,3, 3,3,3, 3,2,1)
  dc_res = dc(X, Y)
  assert_that(are_equal(dc_res$x.to.y, 0))
  assert_that(are_equal(dc_res$y.to.x, 0))
}