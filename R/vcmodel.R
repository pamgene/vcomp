modelOperator = function(df, model, nomfac, reml = TRUE){

  if (length(nomfac) > 0){
    for (i in 1:length(nomfac)){
      df[[ nomfac[i] ]] = as.factor(df[[ nomfac[i] ]])
    }
  }
  models = df %>% group_by(ID, rowSeq) %>% do({
    tryCatch({
      aLme = lmer(model, data = ., REML = reml)
      cdf = data.frame(colSeq = .$colSeq, cValue = fixef(aLme)[1] + resid(aLme), residuals = resid(aLme))
      return(data.table(list(aLme), list(cdf) ))
    },
    error = function(e){
      print(e)
      return(data.table())
    }
    )
  })
}

getVarComp = function(df){
  dfvc = df %>% group_by(ID, rowSeq) %>% do({
    aLme = .$V1[[1]]
    if (!is.null(aLme)){
      vc <- VarCorr(aLme)
      v = as.numeric(vc)
      names(v) = attr(vc, 'names')
      v = c(v, residual = attr(vc, 'sc')^2)
      y0 = fixef(aLme)[1]
      comps = names(v)
      s = sqrt(v)
      cv = sqrt(exp( (s*log(2))^2 ) - 1)
      r = v/sum(v)
      out = data.frame(comp = comps, s = s, r = r, cv = cv, y0 = y0)
    } else {
      out = data.table()
    }
  })
}

getFxdComp = function(df){
  dffxd = df %>% group_by(ID, rowSeq) %>% do({
    aLme = .$V1[[1]]
    if (!is.null(aLme)){
      fxd = fixef(aLme)
      out = data.frame(comp = names(fxd), fxd = fxd, y0 = fxd[1])
    } else {
      out = data.table()
    }
  })
}

getdfout = function(df){
  df_out = df%>%group_by(ID, rowSeq) %>% do({return(.$V2[[1]])})
}
