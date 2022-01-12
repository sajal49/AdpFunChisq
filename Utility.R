# Implementing some utility functions

# Convert a list to matrix (serially)
ListToDataMatrix = function(data_list){
  
  new_data = matrix(0,ncol=length(data_list[[1]]),nrow=length(data_list))
  for(i in 1:nrow(new_data))
    new_data[i,] = data_list[[i]]
  return(new_data)
  
}

# Convert a list to matrix (non-serially)
ListToDataMatrixMethods = function(results){
  
  new_data = matrix(0, nrow=2*length(results),
                    ncol=ncol(results[[1]]))
  for(i in 1:(nrow(new_data)/2))
    new_data[i,] = results[[i]][1,]
  for(i in ((nrow(new_data)/2)+1):(nrow(new_data)))
    new_data[i,] = results[[i-(nrow(new_data)/2)]][2,]
  return(new_data)
  
}

