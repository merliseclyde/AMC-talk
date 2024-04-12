# ASI_sampler function
# converted from matlab by chatgbt
ASI_sampler <- function(data_fixed, data, target, g, gprior, mode, RB_adap, RB, 
                        hparam, tau, nu, burnin, numbofits, thin, numbofreps, 
                        heat, adap_type, updateg) {
  
  # Suppress warnings
  options(warn=-1)
  
  # Initialize variables
  rhostar <- numeric(length(heat) - 1)
  order <- matrix(1:(numbofreps * length(heat)), nrow=numbofreps, ncol=length(heat), byrow=TRUE)
  
  # Calculate rhostar
  for (i in 1:(length(heat) - 1)) {
    rhostar[i] <- log(log(heat[i] / heat[i + 1]))
  }
  
  # Get dimensions of data
  n <- nrow(data)
  p <- ncol(data)
  pfixed <- ncol(data_fixed)
  
  numbofchains <- length(heat)
  
  sigmasqalpha <- 0
  sigmasqbeta <- 0
  
  # Convert single hparam to array if needed
  if (length(hparam) == 1) {
    fixed <- 1
    w <- matrix(hparam, nrow=numbofreps, ncol=numbofchains)
    wstar <- hparam
  } else {
    fixed <- 0
    wa <- hparam[1]
    wb <- hparam[2]
    wstar <- wa / (wa + wb)
  }
  
  logita <- 0.1 / p
  logitb <- 1 - 0.1 / p
  
  XTy_all <- t(data) %*% target
  if (pfixed > 0) {
    XTy_fixed <- t(data_fixed) %*% target
  }
  
  loglike <- matrix(0, nrow=numbofreps, ncol=numbofchains)
  C <- vector("list", length=numbofreps * numbofchains)
  gamma <- array(0, dim=c(numbofreps, p, numbofchains))
  
  for (chain in 1:numbofchains) {
    for (rep in 1:numbofreps) {
      check <- 0
      while (check == 0) {
        gamma[rep, , chain] <- runif(p) < wstar
        
        datastar <- cbind(1, data_fixed, data[, gamma[rep, , chain] == 1])
        
        if (gprior == 1) {
          n0star <- matrix(0, nrow=(1 + size(data_fixed, 2)), ncol=(1 + size(data_fixed, 2)))
          n0star[(1 + size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma[rep, , chain])), (1 + size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma[rep, , chain]))] <- 1 / g[rep, chain] * t(datastar[, 2:end]) %*% datastar[, 2:end]
        } else {
          n0star <- 1 / g[rep, chain] * diag(sum(gamma[rep, , chain]) + 1 + size(data_fixed, 2))
          n0star[1:(1 + size(data_fixed, 2)), 1:(1 + size(data_fixed, 2))] <- 0
        }
        
        C[[rep + (chain - 1) * numbofreps]] <- solve(t(datastar) %*% datastar + n0star)
        
        if (pfixed > 0) {
          XTy <- c(sum(target), XTy_fixed, XTy_all[gamma[rep, , chain] == 1])
        } else {
          XTy <- c(sum(target), XTy_all[gamma[rep, , chain] == 1])
        }
        
        loglike[rep, chain] <- -0.5 * sum(gamma[rep, , chain]) * log(g[rep, chain])
        loglike[rep, chain] <- loglike[rep, chain] - 0.5 * log(det(t(datastar) %*% datastar + n0star))
        loglike[rep, chain] <- loglike[rep, chain] - (0.5 * n + sigmasqalpha) * log(sigmasqbeta + 0.5 * (sum(target^2) - XTy %*% C[[rep + (chain - 1) * numbofreps]] %*% XTy))
        
        if (is.na(loglike[rep, chain]) == FALSE && is.infinite(loglike[rep, chain]) == FALSE) {
          check <- 1
        }
      }
    }
  }
  
  if (fixed == 1) {
    gibbssteps <- c(0, 0)
  } else {
    gibbssteps <- c(0, 0, 0)
  }
  
  numbofstep <- rep(0, numbofchains)
  
  for (rep in 1:numbofreps) {
    for (chain in 1:numbofchains) {
      numbofstep[chain] <- numbofstep[chain] + 1
      
      if (numbofstep[chain] > burnin && numbofstep[chain] %% thin == 0) {
        if (fixed == 1) {
          gibbssteps <- cbind(gibbssteps, w[rep, chain])
        } else {
          gibbssteps <- cbind(gibbssteps, wa / (wa + wb))
        }
      }
      
      if (fixed == 1) {
        w[rep, chain] <- wstar
      } else {
        gamma2 <- gamma[rep, , chain]
        
        XTy <- XTy_all[gamma2 == 1]
        
        if (pfixed > 0) {
          XTy <- c(sum(target), XTy_fixed, XTy)
        } else {
          XTy <- c(sum(target), XTy)
        }
        
        datastar <- cbind(1, data_fixed, data[, gamma2 == 1])
        
        n0star <- 1 / g[rep, chain] * diag(sum(gamma2) + 1 + pfixed)
        n0star[1:(1 + pfixed), 1:(1 + pfixed)] <- 0
        
        C[[rep + (chain - 1) * numbofreps]] <- solve(t(datastar) %*% datastar + n0star)
        
        datastar <- cbind(1, data_fixed, data)
        n0star <- 1 / g[rep, chain] * diag(n + 1 + pfixed)
        n0star[1:(1 + pfixed), 1:(1 + pfixed)] <- 0
        
        Ctot <- solve(t(datastar) %*% datastar + n0star)
        
        gamma3 <- gamma2
        
        for (s in 1:numbofits) {
          if (pfixed > 0) {
            if (updateg == 1) {
              g[rep, chain] <- rgamma(1, nu + 0.5 * (n + sum(gamma2)), tau + 0.5 * (target %*% target - 2 * XTy %*% C[[rep + (chain - 1) * numbofreps]] %*% XTy + trace(Ctot %*% t(datastar) %*% datastar)))
            }
          }
          
          XTy <- XTy_all[gamma3 == 1]
          if (pfixed > 0) {
            XTy <- c(sum(target), XTy_fixed, XTy)
          } else {
            XTy <- c(sum(target), XTy)
          }
          
          datastar <- cbind(1, data_fixed, data[, gamma3 == 1])
          n0star <- 1 / g[rep, chain] * diag(sum(gamma3) + 1 + pfixed)
          n0star[1:(1 + pfixed), 1:(1 + pfixed)] <- 0
          
          C[[rep + (chain - 1) * numbofreps]] <- solve(t(datastar) %*% datastar + n0star)
          
          pstar <- apply(C[[rep + (chain - 1) * numbofreps]], 1, function(x) {
            x %*% XTy
          })
          
          mean1 <- pstar[1]
          means <- pstar[-1]
          
          mean1 <- mean1 / sqrt(sigmasqbeta + 0.5 * (target %*% target - 2 * XTy %*% C[[rep + (chain - 1) * numbofreps]] %*% XTy + trace(Ctot %*% t(datastar) %*% datastar)))
          means <- means / sqrt(sigmasqbeta + 0.5 * (target %*% target - 2 * XTy %*% C[[rep + (chain - 1) * numbofreps]] %*% XTy + trace(Ctot %*% t(datastar) %*% datastar)))
          
          mean1 <- mean1 / sqrt(sigmasqalpha)
          means <- means / sqrt(sigmasqalpha)
          
          u1 <- rnorm(1)
          us <- rnorm(p - pfixed)
          
          temp1 <- qlogis(plogis(mean1 + u1))
          temps <- qlogis(plogis(means + us))
          
          gamma3[gamma3 == 1] <- temps > logitb
          gamma3[gamma3 == 0] <- temps < logita
          
          gamma3[1] <- temp1 > logitb
          gamma3[1] <- temp1 < logita
        }
        
        gamma[rep, , chain] <- gamma3
        
        gibbssteps <- cbind(gibbssteps, wstar)
      }
    }
  }
  
  options(warn=0)
  
  return(gibbssteps)
}
