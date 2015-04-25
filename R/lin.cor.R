lin.corr <- function(x,y, method = "pearson"){
        
        cor.coef <- vector("numeric",1)
                
        if(method == "pearson"){
                cor.coef <-  pearson_corr_centred(x,y)
        }
        
        if(method == "pearson2"){
                cor.coef <-  pearson_corr_noncentred(x,y)
        }
        
        if(method == "sq_pearson"){
                cor.coef <-  squared_pearson_corr(x,y)
        }
        
        return (cor.coef)
}





