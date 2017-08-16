#' @import pgMulticore
#' @export
modelOperator = function(df, model, nomfac, reml = TRUE){

  if (length(nomfac) > 0){
    for (i in 1:length(nomfac)){
      df[[ nomfac[i] ]] = as.factor(df[[ nomfac[i] ]])
    }
  }
  #models = df %>% group_by(rowSeq) %>% do(modelFun(., model, reml))
  aResult = doMultiCore(~rowSeq, df,
                        operatorFunction = modelFun,
                        .packages = c("lme4", "data.table"),
                        model = model,
                        reml = reml)
}

modelFun = function(df, model, reml){

  aLme = try(lmer(model, data = df, REML = reml))

  if(!inherits(aLme, "try-error")){
    aLme = lmer(model, data = df, REML = reml)
    cdf = data.frame(colSeq = df$colSeq, cValue = fixef(aLme)[1] + resid(aLme), residuals = resid(aLme))
    out = data.table(rowSeq = df$rowSeq[1], ID = df$ID[1], aLme = list(aLme), cdf = list(cdf ))
  } else {
    out = data.table(rowSeq = df$rowSeq[1], ID = df$ID[1], aLme = list(), cdf = list())
  }
  return(out)
}

getVarComp = function(df){
  dfvc = df %>% group_by(rowSeq) %>% do({
    aLme = .$aLme[[1]]
    if (!is.null(aLme)){
      vc <- VarCorr(aLme)
      df = adply(names(vc),.margins = 1,  .fun = function(x) {
        comp = vc[[x]]
        comp.sd = attr(comp, "stddev")
        comp.names = paste(x, names(comp.sd), sep = ".")
        result = data.frame(comp = comp.names, s = comp.sd)
      })
      df = subset(df, select = c(comp,s))
      df = rbind(df, data.frame(comp = "residual", s = attr(vc, 'sc')) )
      y0 = fixef(aLme)[1]
      cv = sqrt(exp( (df$s*log(2))^2 ) - 1)
      r = df$s^2/sum(df$s^2)
      df = data.frame(df, r = r, cv = cv, y0 = y0)
    } else {
      df = data.frame()
    }
  })
}

getFxdComp = function(df){
  dffxd = df %>% group_by(rowSeq) %>% do({
    aLme = .$aLme[[1]]
    if (!is.null(aLme)){
      fxd = fixef(aLme)
      data.frame(comp = names(fxd), fxd = fxd, y0 = fxd[1])
    } else {
      data.table()
    }
  })
}

getdfout = function(df){
  df_out = df%>%group_by(ID, rowSeq) %>% do({(.$cdf[[1]])})
}
