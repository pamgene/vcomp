modelOperator = function(df, model){
  model.fac = attr(terms(model), "term.labels")
  if (length(model.fac)>0){
      for (f in colnames(df)){
        if(any(grepl(f, model.fac))) df[[f]] = as.factor(df[[f]])
      }
  }
  models = df %>% group_by(ID, rowSeq) %>% do({
    tryCatch({
      aLme = lmer(model, data = .)
      return(data.table(list(aLme)))
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
      y0 = fixef(aLme)
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


