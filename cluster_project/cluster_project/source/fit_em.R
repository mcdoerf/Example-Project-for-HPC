#####################################################################################
# Julia Wrobel
# February 2025
#
# This file produces implements an EM algorithm for a censored exponential data
#####################################################################################


# lambda0: initial value for lambda
# delta: survival status
# y: survival time or censoring time
em_exp = function(lambda0, delta, y, tol = 1e-9, max_iter = 1000, min_iter = 3){

  lambda_cur = lambda0
  n = length(y)
  iter = 1
  tol_criteria = Inf
  # define vectors to store elements of interest
  # in this case I will store both Q and lambda
  Q = lambda_vec = rep(NA, length = max_iter)

  while(iter < max_iter  & tol_criteria > tol){

    ###############################################################
    ## M-step
    ###############################################################

    lambda_vec[iter] = lambda_cur
    lambda_next = sum(delta*y + (1-delta)*(y + lambda_cur))/n

    ###############################################################

    ###############################################################
    ## E-step
    ###############################################################

    # calculate Q

    Q[iter] = -n * log(lambda_next) - 1/lambda_next * sum(delta*y + (1-delta)*(y + lambda_cur))


    ## check convergence criteria and increase iteration

    if(iter > min_iter){
      tol_criteria = abs(Q[iter] - Q[iter-1])
    }
    iter = iter + 1
    lambda_cur = lambda_next
    #message(paste0("iteration: ", iter, "; ll: ", round(tol_criteria, 4)))
  }

  # final parameter values
  Q[iter] = -n * log(lambda_cur) - 1/lambda_cur * sum(delta*y + (1-delta)*(y+lambda_cur))
  lambda_vec[iter] = lambda_cur

  ### return parameters of interest
  list(solution = lambda_cur,
       iterations = iter,
       Q = Q[1:iter],
       lambda_vec = lambda_vec[1:iter])

}

