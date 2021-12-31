# The original codes for DC was obtained from:
# http://eda.mmci.uni-saarland.de/prj/cisc/

# for testing
require(assertthat)

Cisc = function(x,y)
{
  scX = Stochastic_complexity(x)
  scY = Stochastic_complexity(y)
  
  xgy = makelist(y,x)
  ygx = makelist(x,y)
  
  mXgY = lapply(xgy,function(z)
  {
    return(Stochastic_complexity(z))  
  })
  
  mXgY = sum(unlist(mXgY,use.names = FALSE))
  
  mYgX = lapply(ygx,function(z)
  {
    return(Stochastic_complexity(z))
  })
  
  mYgX = sum(unlist(mYgX, use.names = FALSE))
  
  Cisc.x.to.y = scX + mYgX
  Cisc.y.to.x = scY + mXgY
  
  list(x.to.y = Cisc.x.to.y, y.to.x = Cisc.y.to.x)
}

Stochastic_complexity = function(x)
{
  freqx = table(x)
  n = length(x)
  L = length(freqx)

  logLikely = 0
  
  for(i in freqx)
    logLikely = logLikely + (i*(log2(n) - log2(i)))
  
  nterm = multinomial_with_recurrence(L,n)
  sc = logLikely + log2(nterm) 
  
  return(sc)
}

multinomial_with_recurrence = function(L,n)
{
  total = 1
  b = 1
  d = 10
  
  bound = (ceiling(2 + sqrt(2*n*d*log(10))))

  for(k in 1:bound)
  {
    b = ((n-k+1)*b)/(n)
    total = total + b
  }
  
  old.sum = 1
  if(L>=3)
  {
    for(j in 3:L)
    {
      new.sum = total + ((n*old.sum)/(j-2))
      old.sum = total
      total = new.sum
    } 
  }
  
  if(L==1)
    total = 1
  
  return(total)
}

makelist = function(x,y)
{
  unx = unique(x)
  ylist = list()
  
  for(i in 1:length(unx))
  {
    ylist[[i]] = y[x==unx[i]]
  }
  
  return(ylist)
}


test_cisc = function(){
  
 X = c(1,1,1,1, 2,2,2,2, 3,3,3,3)
 Y = c(1,1,1,1, 2,2,2,2, 1,1,1,1)
 res = Cisc(X, Y)
 assert_that(are_equal(res$x.to.y, 23.110070977710365))
 assert_that(are_equal(res$y.to.x, 23.437620070757482))
 
 X = c(1,1,1, 2,2,2, 3,3,3, 1,1,1, 4,4,4,4)
 Y = c(2,2,2, 1,1,1, 2,2,2, 2,2,2, 1,1,1,1)
 res = Cisc(X, Y)
 assert_that(are_equal(res$x.to.y, 37.0538059050482))
 assert_that(are_equal(res$y.to.x, 37.65551544822503))
 
 X = c(1,2,3, 3,2,1, 1,1,1, 2,2,2, 3,3,3)
 Y = c(1,1,1, 2,2,2, 1,2,3, 3,3,3, 3,2,1)
 res = Cisc(X, Y)
 assert_that(are_equal(res$x.to.y, 59.47692437100562))
 assert_that(are_equal(res$y.to.x, 59.47692437100562))
 
 X = c(1,2,3, 3,2,1, 1,1,1, 2,2,2, 3,3,3)
 Y = c(3,2,1, 1,2,3, 3,3,3, 1,1,1, 2,2,2)
 res = Cisc(X, Y)
 assert_that(are_equal(res$x.to.y, 41.46798642990164))
 assert_that(are_equal(res$y.to.x, 41.46798642990164))
 
 X = c(1, 2, 3, 3, 2, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3)
 Y = c(2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1)
 res = Cisc(X, Y)
 assert_that(are_equal(res$x.to.y, 44.39903733951914))
 assert_that(are_equal(res$y.to.x, 44.114233535914465))
 
}
