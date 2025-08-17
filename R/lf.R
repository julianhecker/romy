

#' ROMY-LF
#'
#' This function performs the conditional latent factor model computations
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param epsilon Convergence threshold, default value is 0.01.
#' @param K_min Minimum number of latent factors, default value is 1.
#' @param K_max Maximum number of latent factors, default value is 5.
#' @param max_iter Maximum number of iterations, default value is 20.
#' @param eta Eta parameter, default is value is 0.5.
#'
#' @export
romy_lf=function(Y, X, inds_X, epsilon=0.01, K_min=1, K_max=5, max_iter=20, eta=0.5)
{
  
  N=nrow(Y)
  M=ncol(Y)
  results=list()
  Z_hat_truncated=NULL
  Z_hat_untruncated=NULL
  
  num=20
  r <- (10 / 0.1)^(1/(num-1))
  seqs_scaling_dec <- 1/(0.1 * r^ (0:(num-1)))
  seqs_scaling_inc <- 0.1 * r^ (0:(num-1))
  
  convergence=FALSE
  degenerated=TRUE
  
  for(c_lambda in seqs_scaling_dec)
  {
		res=.romy_lf_algorithm_gamma(Y=Y, X=X, inds_X=inds_X, epsilon=epsilon, K_max=K_max, 
                          c_lambda=c_lambda, max_iter=max_iter, eta=eta)
		
		if(res$convergence & !res$all_zero)
		{
			K_ests=rep(0, length(seqs_scaling_inc))
			ctr=1
			for(c_delta in seqs_scaling_inc)
			{
				K_est=.estimate_K(res$ev_values, c_delta, N=N, M=M)
				if(K_est>=K_min & K_est<=K_max)
				{
					Z_hat_untruncated=res$Z_hat[,1:K_est]
				}
				K_ests[ctr]=K_est; ctr=ctr+1
			}
			if(max(K_ests)<=K_min){
				Z_hat_truncated=res$Z_hat[,1:K_min]
			}
			if(min(K_ests)>=K_max){
				Z_hat_truncated=res$Z_hat[,1:K_max]
			}
			
			convergence=TRUE
			degenerated=FALSE
		}
  }
  
  Z_hat=as.matrix(Z_hat_truncated)
  if(!is.null(Z_hat_untruncated))
  {
    Z_hat=as.matrix(Z_hat_untruncated)
  }
  K_hat=0
  if(!is.null(Z_hat))
  {
	K_hat=ncol(Z_hat)
  }
  
  return(list(Z_hat=Z_hat, K_hat=K_hat, convergence=convergence, degenerated=degenerated))
}

#' Estimates K as described by Chen (2022).
#'
#' @param ev_values Eigenvalues from projected Gamma Matrix product
#' @param c_delta Constant for Delta Computation
#' @param N Sample Size
#' @param M Number of Phenotypes
#'
.estimate_K=function(ev_values, c_delta, N, M)
{
	delta_NM=c_delta*(N+M)*log(N)
	inds=1:N
	K_hat=0
	if(sum(ev_values>=delta_NM)>=1)
	{
		  K_hat=max(inds[ev_values>=delta_NM])
	}
	return(K_hat)
}


#' Performs the proximal gradient descent algorithm as described by Chen 2022 (Appendix B).
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param epsilon Convergence threshold, default value is 0.01.
#' @param c_lambda factor for Lambda.
#' @param max_iter Maximum number of iterations, default value is 20.
#' @param eta Eta parameter, default is value is 0.5.
#' @param K_max Maximum number of latent factors, default value is 5.
#'
.romy_lf_algorithm_gamma=function(Y, X, inds_X, epsilon, c_lambda, max_iter=20, eta=0.5, K_max=5)
{
  N=nrow(Y)
  M=ncol(Y)
  P=sum(lengths(inds_X))
  
  
  w_km1=1
  w_k=1
  
  Gamma_s_km1=matrix(rnorm(N*P), nrow=N, ncol=P)
  Gamma_s_k=matrix(rnorm(N*P), nrow=N, ncol=P)
  
  #################
  L=-1
  for(m in 1:M)
  {
    tmp_X=as.matrix(X[,inds_X[[m]]])
    tmp_X_sqsum=rowSums(tmp_X**2)
    if(L<max(tmp_X_sqsum))
    {
      L=max(tmp_X_sqsum)
    }
  }
  #################
  
  lambda_NM=c_lambda*sqrt((N+M)*log(N))
  
  
  
  tau_0=L
  tau_km1=tau_0
  
  
  k=1
  
  convergence=FALSE
  
  while(!convergence & k<=max_iter){
		# step 1
		
		Gamma_k = Gamma_s_k + (w_km1-1)/w_k *(Gamma_s_k - Gamma_s_km1)
		
		# step 2
		tau_hat_j = eta * tau_km1
		run=TRUE
		iter_in=0
		while(run & iter_in<=max_iter)
		{
		  iter_in=iter_in+1
		  #print(tau_hat_j)
		  x = 1/tau_hat_j * lambda_NM
		  A = Gamma_k - 1/tau_hat_j*.gradient_f_func(Y, X, inds_X, Gamma_k)
		  A_j = .S_operator(A, x)
		  if(is.null(A_j)){Z_hat=0; Gamma_hat=0; ev_values=0; all_zero=FALSE; return(list(convergence=convergence, Z_hat=as.matrix(Z_hat), Gamma_hat=as.matrix(Gamma_hat), ev_values=ev_values, all_zero=all_zero))}
		  Ftmp=.cF_func(Y, X, inds_X, A_j)
		  if(is.null(Ftmp)){Z_hat=0; Gamma_hat=0; ev_values=0; all_zero=FALSE; return(list(convergence=convergence, Z_hat=as.matrix(Z_hat), Gamma_hat=as.matrix(Gamma_hat), ev_values=ev_values, all_zero=all_zero))}
		  Qtmp=.Q_func(Y, X, inds_X, A_j, Gamma_k, tau_hat_j, lambda_NM)
		  if(is.null(Qtmp)){Z_hat=0; Gamma_hat=0; ev_values=0; all_zero=FALSE; return(list(convergence=convergence, Z_hat=as.matrix(Z_hat), Gamma_hat=as.matrix(Gamma_hat), ev_values=ev_values, all_zero=all_zero))}
		  if(Ftmp<=Qtmp)
		  {
			tau_k = tau_hat_j
			run = FALSE
			break
		  }
		  else{
			tau_hat_j = min(1/eta*tau_hat_j, tau_0)
		  }
		}
		tau_k=tau_0
		# step 3
		x = 1/tau_k * lambda_NM
		A = Gamma_k - 1/tau_k * .gradient_f_func(Y, X, inds_X, Gamma_k)
		
		Gamma_s_kp1 = .S_operator(A, x)
		if(is.null(Gamma_s_kp1)){Z_hat=0; Gamma_hat=0; ev_values=0; all_zero=FALSE; return(list(convergence=convergence, Z_hat=as.matrix(Z_hat), Gamma_hat=as.matrix(Gamma_hat), ev_values=ev_values, all_zero=all_zero))}
		
		
		# step 4
		w_kp1 = (1 + sqrt(1 + 4*w_k**2))/2
		
		# step 5
		
		D_kp1= tau_k *(Gamma_k - Gamma_s_kp1) + .gradient_f_func(Y, X, inds_X, Gamma_s_kp1) - .gradient_f_func(Y, X, inds_X, Gamma_k)
		norm_D_kp1=sum(D_kp1**2)
		factor_value=tau_k*max(1, sum(Gamma_s_kp1**2))
		#cat("norm/factor:", norm_D_kp1/factor, "\n")
		if(norm_D_kp1/factor_value<=epsilon)
		{
		  Gamma_out=Gamma_s_kp1
		  convergence=TRUE
		  break
		}
		else{
		  w_km1=w_k
		  w_k=w_kp1
		  tau_km1=tau_k
		  Gamma_s_km1=Gamma_s_k
		  Gamma_s_k=Gamma_s_kp1
		  k=k+1
		}
    
  }
  
  if(convergence)
  {
		Gamma_hat=Gamma_out
		if(var(as.vector(Gamma_hat))>10^-6){
			M_N=(diag(N)-1/N*matrix(1, nrow=N, ncol=N))
			M_Gamma=M_N %*%Gamma_hat
			M_GammaGamma_M=M_Gamma %*% t(M_Gamma)
			evd=eigen(M_GammaGamma_M)
			inds=1:N
			Z_hat=evd$vectors[,1:K_max]*sqrt(N)
			ev_values=evd$values
			all_zero=FALSE
		}else{
			Z_hat=0
			ev_values=0
			all_zero=TRUE
		}
  }else{
		Z_hat=0
		Gamma_hat=0
		ev_values=0
		all_zero=FALSE
  }
  
  return(list(convergence=convergence, Z_hat=as.matrix(Z_hat), Gamma_hat=as.matrix(Gamma_hat), ev_values=ev_values, all_zero=all_zero))
}





#
#' Computes nuclear norm of Matrix Gamma
#'
#' @param Gamma Matrix
#'
.nuclear_norm <- function(Gamma) 
{
  svd_tmp <- tryCatch(
    svd(Gamma),
    error = function(e) {
      message("SVD failed: ", e$message)
      return(NULL)  
    }
  )
  
  if (is.null(svd_tmp)) {
    return(NULL)
  }
  
  return(sum(svd_tmp$d))
}

#' f function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#'
.f_func=function(Y, X, inds_X, Gamma)
{
	M=ncol(Y)
    Yp=matrix(0, nrow=nrow(Y), ncol=ncol(Y))
    for(m in 1:M)
    {
      Yp[,m]=rowSums(as.matrix(X[,inds_X[[m]]]) * as.matrix(Gamma[,inds_X[[m]]]))
    }
    diff=Y-Yp
    f=sum(diff**2)
    return(f)
}

#' gradient f function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#'
.gradient_f_func=function(Y, X, inds_X, Gamma)
{
    N=nrow(Y)
	M=ncol(Y)
    P=ncol(Gamma)
    grad=matrix(0, nrow=N, ncol=P)
    for(m in 1:M)
    {
      pred=rowSums(as.matrix(X[,inds_X[[m]]]) * as.matrix(Gamma[,inds_X[[m]]]))
      grad[, inds_X[[m]]]=-(Y[,m]-pred)*X[,inds_X[[m]]]
    }
    return(grad)
}

#' F function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#'
.cF_func=function(Y, X, inds_X, Gamma)
{
    f=.f_func(Y, X, inds_X, Gamma)
    nn=.nuclear_norm(Gamma)
	if(!is.null(nn)){
		return(f+nn)
	}
	else{
		return(NULL)
	}
}

#' Q function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#' @param Gamma_k Matrix Gamma_k
#' @param tau_k Parameter tau_k
#' @param lambda Parameter lambda
#'
.Q_func=function(Y, X, inds_X, Gamma, Gamma_k, tau_k, lambda)
{
    f=.f_func(Y, X, inds_X, Gamma_k)
    grad=.gradient_f_func(Y, X, inds_X, Gamma_k)
    trace= sum((Gamma-Gamma_k) * grad)
    diff=tau_k/2 * sum((Gamma-Gamma_k)**2)
	nn=.nuclear_norm(Gamma)
	if(!is.null(nn)){
		penalty=lambda*nn
		return(f+trace+diff+penalty)
	}
	else{
		return(NULL)
	}
}

#' S operator as described in Appendix B of Chen 2022.
#'
#' @param A Matrix A.
#' @param x Threshold x.
#'
.S_operator=function(A, x)
{
    
	svd_tmp <- tryCatch(
		svd(A),
		error = function(e) {
		  message("SVD failed: ", e$message)
		  return(NULL)  
		}
	)
	  
	if (is.null(svd_tmp)) {
		return(NULL)
	}else{
		diag_mat=diag(ifelse(svd_tmp$d-x>0, svd_tmp$d-x,0))
		A_new=svd_tmp$u %*% diag_mat %*% t(svd_tmp$v)
		return(A_new)
	}
}



