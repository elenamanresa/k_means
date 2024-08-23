require(data.table)

# Simulate group data and estimate it
# Functions:
# (1) generate data
# (2) Lloyds algorithm

gen_data <-function(dgp,prob_G){
  Ndim = dgp@Ndim
  Tdim = dgp@Tdim
  G = dgp@G
  siga = dgp@sig_alpha
  sigerr = dgp@sig_err
  
  probabilities <- seq(1/(G+1),1,1/G)
  alpha <-  qnorm(probabilities, mean = 0, sd = siga)
  error <- matrix(rnorm(Ndim*Tdim,0,sigerr),Ndim,Tdim)
  assign <- sample(sample(1:G,Ndim,replace=TRUE,prob=prob_G))
  alpha_mat <- matrix(rep(alpha[assign],Tdim),Ndim,Tdim)
  y = alpha_mat + error
  return(list(y,alpha,assign))
}

# (2) Lloyd's algorithm
Lloyds <- function(dgp,y,alpha_ini,itermax){
  flag <- 0
  Ndim <- dgp@Ndim
  Tdim <- dgp@Tdim
  G <- dgp@G
  tol <- 0.0001
  dif <- 1
  data <- CJ(i=1:Ndim,t=1:Tdim)
  data[,y:=as.vector(t(y))]
  
  res = matrix(0,Ndim,G)
  alpha_new = rep(0,G)
  iter = 1
  while(dif > tol &&  iter < 50){
    # assignment
    for(gg in 1:G){
      resaux = y -  matrix(rep(alpha_ini[gg],Ndim*Tdim),Ndim,Tdim)
      resaux2 = resaux^2
      res[,gg] = rowSums(resaux2)
    }
    new_assign = apply(res, 1, FUN = which.min)
    # check if all groups are full
    if(length(unique(new_assign))<G){
      iter = 50
      flag = -1
      #print(iter)
    }else{
    data[,assign:=as.character(new_assign[i])]
    # update alpha
    model <- glm(formula = y ~ assign + 0, data = data)
    dif = 0
    for(gg in 1:G){
      alpha_new[gg] = model$coefficients[[gg]]
      dif = dif +  (alpha_new[gg] - alpha_ini[gg])^2
    }
    #print(dif)
    alpha_ini = alpha_new
    iter = iter+1
    #print(iter)
    data[,alpha_out:=alpha_ini[as.numeric(assign)]]
    }
  }
  return(list(data,flag))
}


# MAIN

setClass("DGP", slots = list(Ndim = "numeric",Tdim = "numeric",G = "numeric",sig_alpha = "numeric", sig_err = "numeric"))


dgp <- new("DGP",Ndim=1000,Tdim=50,G=5,sig_alpha = 1,sig_err=1)
prob_G = rep(1/dgp@G,dgp@G)
#prob_G = c(0.5,0.5)
itermax = 50

# Generate data:
data_out <- gen_data(dgp,prob_G)
y <- data_out[[1]]
alpha0 <- data_out[[2]]
assign0 <- data_out[[3]]

# Estimate the model
inicond = 100
min_loss_mat = matrix(0,inicond,1)
for(ii in 1:inicond){
  print(ii)
  alpha_ini = sort(quantile(y,seq(1/(dgp@G+1),1-1/(dgp@G+1),1/(dgp@G+1))) + rnorm(dgp@G))
  LLoydsout <- Lloyds(dgp,y,alpha_ini)
  dataoutput <- LLoydsout[[1]]
  flag <- LLoydsout[[2]]
  if (flag == 0){
    loss = sum((dataoutput$y-dataoutput$alpha_out)^2)/(dgp@Ndim*dgp@Tdim)
  }else{
    loss =dgp@sig_err*dgp@Ndim*dgp@Tdim
  }
  missclass <- 1-sum(dataoutput$assign[seq(1, dgp@Ndim*dgp@Tdim, by = dgp@Tdim)] == assign0)/dgp@Ndim
  print(list(missclass,loss))
  if(ii ==1){
    min_loss = loss
    best_data = dataoutput
    #missclass <- 1-sum(best_data$assign[seq(1, dgp@Ndim*dgp@Tdim, by = dgp@Tdim)] == assign0)/dgp@Ndim
    #print(list(missclass,loss))
    best = ii
  }else{ 
    if(loss < min_loss){
      min_loss = loss
      best_data = dataoutput
      # Compute missclassification probability:
      best = ii
    }
  }
  min_loss_mat[ii] = min_loss

}



missclass <- 1-sum(best_data$assign[seq(1, dgp@Ndim*dgp@Tdim, by = dgp@Tdim)] == assign0)/dgp@Ndim
print(list(best,missclass,min_loss))

if(inicond > 1){
# plot loss
plot(min_loss_mat[min_loss_mat<1000])
}
# plot alpha_hat vs alpha0
print("alpha_hat")
print(c(sort(unique(best_data$alpha_out))))
print("alpha true")
alpha0


