## The file contains R functions that deals with data from Immunohistochemistry staining. This type of data frame contains staining score and various patient pathological information. The purpose of this anaylsis is to review any association between tissue staining score with the information.


## The following function performs this by pulling "bin" (dichotomized staining score) with "ind" (patient indicators).
## ordinal_kendall is a vector with each element corresponding to column number of "ind" that undergo kendall test. 
## recur_cox is a vector of length 2;  first element is recurrence time and second element is the censor.
## mort_cox is a vector of length 3; first and 2nd elements are the mortality time and survival after recurrence. 3rd element is the censor. 


stain_Func <- function(bin, ind, ordinal_kendall, recur_cox, mort_cox, bin_fisher = setdiff(seq(1:length(indicator)), c(ordinal_kendall, recur_cox, mort_cox))){
    pointer_mapply <- expand.grid(1:ncol(bin), 1:(ncol(ind)-2))
    mapply(function(b, i){
        if (i %in% ordinal_kendall){
            Kendall(bin[[b]], ind[[i]])$sl
        } else if (i %in% recur_cox[1]) {
            if (length(which(bin[[b]] == TRUE))>=3 & length(which(bin[[b]] == FALSE))>=3){
                coxph(Surv(ind[[i]], ind[[recur_cox[2]]])~bin[[b]]) %>% summary(.) %>% .$coefficients %>% .[[5]]
            } else {return("NA")}
        } else if (i %in% mort_cox[1:2]) {
            if (length(which(bin[[b]] == TRUE))>=3 & length(which(bin[[b]] == FALSE))>=3){
                coxph(Surv(ind[[i]], ind[[mort_cox[3]]])~bin[[b]]) %>% summary(.) %>% .$coefficients %>% .[[5]] # use the pval for staining cutoff, not the likelihood pval for the model fitted 
            } else {return("NA")}
        } else {
            fisher.test(as.factor(bin[[b]]), as.factor(ind[[i]]))$p.value
        }
    }, pointer_mapply[,1], pointer_mapply[,2]) %>% 
        as.numeric(.) %>% 
        matrix(., nrow = ncol(bin), byrow = FALSE) %>% 
        tbl_df()
}

stain_coef_Func <- function(bin, ind, ordinal_kendall, recur_cox, mort_cox, bin_fisher = setdiff(seq(1:length(indicator)), c(ordinal_kendall, recur_cox, mort_cox))){
    pointer_mapply <- expand.grid(1:ncol(bin), 1:(ncol(ind)-2))
    mapply(function(b, i){
        if (i %in% ordinal_kendall){
            Kendall(bin[[b]], ind[[i]])$tau
        } else if (i %in% recur_cox[1]) {
            if (length(which(bin[[b]] == TRUE))>=3 & length(which(bin[[b]] == FALSE))>=3){
                coxph(Surv(ind[[i]], ind[[recur_cox[2]]])~bin[[b]]) %>% summary(.) %>% .$coefficients %>% .[[2]]
            } else {return("NA")}
        } else if (i %in% mort_cox[1:2]) {
            if (length(which(bin[[b]] == TRUE))>=3 & length(which(bin[[b]] == FALSE))>=3){
                coxph(Surv(ind[[i]], ind[[mort_cox[3]]])~bin[[b]]) %>% summary(.) %>% .$coefficients %>% .[[2]] # use the pval for staining cutoff, not the likelihood pval for the model fitted 
            } else {return("NA")}
        } else {
            fisher.test(as.factor(bin[[b]]), as.factor(ind[[i]]))$estimate
        }
    }, pointer_mapply[,1], pointer_mapply[,2]) %>% 
        as.numeric(.) %>% 
        matrix(., nrow = ncol(bin), byrow = FALSE) %>% 
        tbl_df()
}










## ## for xtable 
paste_coef_Func <- function(pval_tbl, coef_tbl){
    combi <- sapply(1:length(as.matrix(pval_tbl)), function(n){
        paste0(signif(as.numeric(as.matrix(pval_tbl)[n]), digits = 3), "\n(", signif(as.numeric(as.matrix(coef_tbl)[n]), digits = 3), ")")
    })
    sapply(1:length(combi), function(n){ifelse(as.numeric(unlist(strsplit(combi[n], "\\n"))[1])<0.1, ifelse(as.numeric(unlist(strsplit(combi[n], "\\n"))[1])<0.05, paste0("\\textbf{", combi[n], "}"), combi[n]), "")}) %>%
        matrix(., nrow = nrow(pval_tbl)) %>% tbl_df()
}


