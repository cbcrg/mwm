predict.censReg <- function(model, data, type = "fitted") { 
     estim <- model$estimate
     coefficients <- estim[-grep("logSigma", names(estim))]
     names.coef <- names(coefficients) 
     indices <- which(names(coefficients)%in%colnames(data))

     if (length(indices)==(length(coefficients))) { 
        fitted <- (data %*% coefficients)[,1]
     } 
     else { 
      stop("The dataset does not contain some of the variables.") 
     } 

     sigma <- summary(model)$estimate[length(summary(model)$estimate[,1]),1] 

     if (type == "fitted") { 
       fitted 
     } 
     else if ((type == "censored") || (type == "truncated")) { 
       if (model$left == -Inf) { 
        print("Right censored/truncated mean:") 
        right <- model$right 
        truncated <- fitted - sigma*(lambda((fitted - right)/sigma)) 
        prob_trunc <- 1 - pnorm((fitted - right)/sigma) 
        censored <- right*dnorm((fitted - right)/sigma) + prob_trunc*truncated 
       } 
       else if (model$right == Inf) { 
        print("Left censored/truncated mean:") 
        left <- model$left 
        truncated <- fitted + sigma*(lambda((fitted - left)/sigma)) 
        prob_trunc <- pnorm((fitted - left)/sigma) 
        censored <- left*dnorm((fitted - left)/sigma) + prob_trunc*truncated 
       } 
       else { 
        print("Interval censored/truncated mean:") 
        left <- model$left 
        right <- model$right
        truncated <- fitted - sigma*((dnorm((left-fitted)/sigma)-dnorm((right-fitted)/sigma))/(pnorm((left-fitted)/sigma)-
                       pnorm((right-fitted)/sigma))) 
        prob_trunc <- pnorm((fitted - right)/sigma) - pnorm((fitted - left)/sigma) 
        censored <- left*dnorm((fitted - left)/sigma) + prob_trunc*truncated + right*dnorm((fitted - right)/sigma)
       } 
       if (type == "censored") { 
        censored 
       } 
       else { 
        truncated 
       }
      } 
} 


lambda <- function(x) { 
   lambda <- dnorm(x)/pnorm(x) 
   lambda 
}

