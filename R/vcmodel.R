modelOperator = function(df, model, nomfac, reml = TRUE){

  if (length(nomfac) > 0){
    for (i in 1:length(nomfac)){
      df[[ nomfac[i] ]] = as.factor(df[[ nomfac[i] ]])
    }
  }
  models = df %>% group_by(ID, rowSeq) %>% do(modelFun(., model, reml))
}

modelFun = function(., model, reml){
  aLme = try(lmer(model, data = ., REML = reml))

  if(!inherits(aLme, "try-error")){
    aLme = lmer(model, data = ., REML = reml)
    cdf = data.frame(colSeq = .$colSeq, cValue = fixef(aLme)[1] + resid(aLme), residuals = resid(aLme))
    out = data.table(list(aLme), list(cdf) )
  } else {
    out = data.table()
  }
  return(out)
}

getVarComp = function(df){
  dfvc = df %>% group_by(ID, rowSeq) %>% do({
    aLme = .$V1[[1]]
    if (!is.null(aLme)){
      vc <- VarCorr(aLme)
      df = adply(names(vc),.margins = 1,  .fun = function(x) {
        comp = vc[[x]]
        comp.sd = attr(comp, "stddev")
        comp.names = paste(x, names(comp.sd), sep = ".")
        result = data.frame(comp = comp.names, s = comp.sd)
      })
      df = subset(df, select = c(comp,s))
      df = rbind(df, data.frame(comp = "residual", s = c(df$s,attr(vc, 'sc')) ))
      y0 = fixef(aLme)[1]
      cv = sqrt(exp( (df$s*log(2))^2 ) - 1)
      r = df$s^2/sum(df$s^2)
      df = data.frame(df, r = r, cv = cv, y0 = y0)
    } else {
      data.table()
    }
  })
}

getFxdComp = function(df){
  dffxd = df %>% group_by(ID, rowSeq) %>% do({
    aLme = .$V1[[1]]
    if (!is.null(aLme)){
      fxd = fixef(aLme)
      data.frame(comp = names(fxd), fxd = fxd, y0 = fxd[1])
    } else {
      data.table()
    }
  })
}

getdfout = function(df){
  df_out = df%>%group_by(ID, rowSeq) %>% do({(.$V2[[1]])})
}
