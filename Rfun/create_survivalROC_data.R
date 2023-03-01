create_survivalROC_data <- 
function (tvec, data, marker, time = time, status = status) 
{
    message("Rfun/create_survivalROC_data(), March 01, 2023")
    dnm <- as.character(substitute(data))
    mm <- as.character(substitute(marker))
    tt <- as.character(substitute(time))
    ss <- as.character(substitute(status))
    cnms <- colnames(data)
    exit <- FALSE
    if (!any(cnms == mm)) 
        exit <- TRUE
    if (!any(cnms == tt)) 
        exit <- TRUE
    if (!any(cnms == ss)) 
        exit <- TRUE
    if (exit) 
        stop("One(or more) variables: ", mm, ", ", 
            tt, ", ", ss, " not found in ", dnm)
    tibble(t = tvec) %>% mutate(survivalROC = map(t, survivalROC_helper, 
        data = data, mm = mm, tt = tt, ss = ss), auc = map_dbl(survivalROC, 
        magrittr::extract2, "AUC"), df_survivalROC = map(survivalROC, 
        function(obj) {
            as_tibble(obj[c("cut.values", "TP", "FP")])
        })) %>% select(-survivalROC) %>% unnest(df_survivalROC) %>% 
        arrange(t, FP, TP)
}
