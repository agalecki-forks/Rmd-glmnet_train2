survivalROC_helper <-
function (t, data, mm, tt, ss) {
    message("survivalROC_helper(), March 1, 2023")
    survivalROC::survivalROC(Stime = data[[tt]], status = data[[ss]], 
        marker = data[[mm]], predict.time = t, method = "NNE", 
        span = 0.25 * nrow(data)^(-0.2))
}
