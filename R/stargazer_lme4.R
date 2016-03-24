
#' stargazer function for lme4 with random effects
#' @export
stargazer_lme4 <- function(..., type = "latex", title = "", style = "default", 
                           summary = NULL, out = NULL, out.header = FALSE, column.labels = NULL, 
                           column.separate = NULL, covariate.labels = NULL, dep.var.caption = NULL, 
                           dep.var.labels = NULL, dep.var.labels.include = TRUE, align = FALSE, 
                           coef = NULL, se = NULL, t = NULL, p = NULL, t.auto = TRUE, 
                           p.auto = TRUE, ci = FALSE, ci.custom = NULL, ci.level = 0.95, 
                           ci.separator = NULL, add.lines = NULL, apply.coef = NULL, 
                           apply.se = NULL, apply.t = NULL, apply.p = NULL, apply.ci = NULL, 
                           colnames = NULL, column.sep.width = "5pt", decimal.mark = NULL, 
                           df = TRUE, digit.separate = NULL, digit.separator = ",", 
                           digits = 2, digits.extra = NULL, flip = FALSE, float = TRUE, 
                           float.env = "table", font.size = NULL, header = TRUE, initial.zero = NULL, 
                           intercept.bottom = TRUE, intercept.top = FALSE, keep = NULL, 
                           keep.stat = NULL, label = "", model.names = NULL, model.numbers = TRUE, 
                           multicolumn = TRUE, no.space = NULL, notes = NULL, notes.align = NULL, 
                           notes.append = TRUE, notes.label = NULL, object.names = FALSE, 
                           omit = NULL, omit.labels = NULL, omit.stat = NULL, omit.summary.stat = NULL, 
                           omit.table.layout = NULL, omit.yes.no = c("Yes", "No"), order = NULL, 
                           ord.intercepts = FALSE, perl = FALSE, report = NULL, rownames = NULL, 
                           rq.se = "nid", selection.equation = FALSE, single.row = FALSE, 
                           star.char = NULL, star.cutoffs = NULL, suppress.errors = FALSE, 
                           table.layout = NULL, table.placement = "!htbp", zero.component = FALSE, 
                           summary.logical = TRUE, summary.stat = NULL, nobs = TRUE, 
                           mean.sd = TRUE, min.max = TRUE, median = FALSE, iqr = FALSE) {
  
  
  Tables <- capture.output(stargazer(..., 
                                     type = type, 
                                     title = title, 
                                     style = style, 
                                     summary = summary, 
                                     out = out, 
                                     out.header = out.header, 
                                     column.labels = column.labels, 
                                     column.separate = column.separate, 
                                     covariate.labels = covariate.labels, 
                                     dep.var.caption = dep.var.caption, 
                                     dep.var.labels = dep.var.labels, 
                                     dep.var.labels.include = dep.var.labels.include, 
                                     align = align, 
                                     coef = coef, 
                                     se = se, 
                                     t = t, 
                                     p = p, 
                                     t.auto = t.auto, 
                                     p.auto = p.auto, 
                                     ci = ci, 
                                     ci.custom = ci.custom, 
                                     ci.level = ci.level, 
                                     ci.separator = ci.separator, 
                                     add.lines = add.lines, 
                                     apply.coef = apply.coef, 
                                     apply.se = apply.se, 
                                     apply.t = apply.t, 
                                     apply.p = apply.p, 
                                     apply.ci = apply.ci, 
                                     colnames = colnames, 
                                     column.sep.width = column.sep.width, 
                                     decimal.mark = decimal.mark, 
                                     df = df, 
                                     digit.separate = digit.separate, 
                                     digit.separator = digit.separator, 
                                     digits = digits, 
                                     digits.extra = digits.extra, 
                                     flip = flip, 
                                     float = float, 
                                     float.env = float.env, 
                                     font.size = font.size, 
                                     header = header, 
                                     initial.zero = initial.zero, 
                                     intercept.bottom = intercept.bottom, 
                                     intercept.top = intercept.top, 
                                     keep = keep, 
                                     keep.stat = keep.stat, 
                                     label = label, 
                                     model.names = model.names, 
                                     model.numbers = model.numbers, 
                                     multicolumn = multicolumn, 
                                     no.space = no.space, 
                                     notes = notes, 
                                     notes.align = notes.align, 
                                     notes.append = notes.append, 
                                     notes.label = notes.label, 
                                     object.names = object.names, 
                                     omit = omit, 
                                     omit.labels = omit.labels, 
                                     omit.stat = omit.stat, 
                                     omit.summary.stat = omit.summary.stat, 
                                     omit.table.layout = omit.table.layout, 
                                     omit.yes.no = omit.yes.no, 
                                     order = order, 
                                     ord.intercepts = ord.intercepts, 
                                     perl = perl, 
                                     report = report, 
                                     rownames = rownames, 
                                     rq.se = rq.se, 
                                     selection.equation = selection.equation, 
                                     single.row = single.row, 
                                     star.char = star.char, 
                                     star.cutoffs = star.cutoffs, 
                                     suppress.errors = suppress.errors, 
                                     table.layout = table.layout, 
                                     table.placement = table.placement, 
                                     zero.component = zero.component, 
                                     summary.logical = summary.logical, 
                                     summary.stat = summary.stat, 
                                     nobs = nobs, 
                                     mean.sd = mean.sd, 
                                     min.max = min.max, 
                                     median = median, 
                                     iqr = iqr)
  )
  
  models <- list(...)
  n.models <- length(models)
  
  Tables <- as.data.frame(Tables)
  Tables$Tables <- as.character(Tables$Tables)
  
  fixed.idx <- grep("Constant", Tables$Tables) - 2L
  
  if (!is.null(column.labels) & !model.numbers) {
    Tables$Tables[fixed.idx] <- gsub(paste("&", paste(column.labels, collapse = " & ")), 
                                     paste("{\\\\bf Fixed Effects} &", paste(column.labels, collapse = " & ")), 
                                     Tables$Tables[fixed.idx])
  } else {
    Tables$Tables[fixed.idx] <- gsub("\\[-1.8ex]", "\\[-1.8ex] {\\\\bf Fixed Effects}", Tables$Tables[fixed.idx])
  }
  
  ## Replace the use of the term 'Constant' with intercept
  Tables$Tables <- gsub("Constant", "(Intercept)", Tables$Tables)
  
  if ("s" %in% omit.table.layout) {
    r <- nrow(Tables) - 5
  } else {
    r <- nrow(Tables) - 10
  }
  randomeffect <- paste0("{\\bf Random Effects} & ", paste(rep("Std.Dev.", n.models), collapse = " & "), "\\\\")
  hline <- "\\hline\\\\[-1.8ex]"
  newline <- "\\\\"
  Tables <- insertrow(Tables, hline, r)
  Tables <- insertrow(Tables,randomeffect, r+1)
  Tables <- insertrow(Tables, hline, r+2)
  
  added.row.counter <- 2
  
  ## get all unique random effect names across all models
  all.ranef.names <- all.grp.names <- NULL
  for (m in 1:n.models) {
    vc <- VarCorr(models[[m]])
    reStdDev <- c(lapply(vc, attr, "stddev"))
    repLens <- lengths(reStdDev)
    grp.names <- names(vc)
    cur.names <- rep(NA, sum(repLens))
    name.names <- unlist(lapply(vc, colnames))
    all.grp.names <- c(grp.names, all.grp.names)
    
    k <- 0
    for (g in 1:length(grp.names)) {
      for (gs in 1:repLens[g]) {
        k <- k + 1
        
        if (gs == 1) {
          cur.names[k] <- paste(grp.names[g], name.names[k], sep = " - ")
        } else {
          cur.names[k] <- paste0(grp.names[g], " - ", name.names[k])
        }
      }
    }
    
    all.ranef.names <- c(all.ranef.names, cur.names)
  }
  all.ranef.names <- sort(unique(all.ranef.names))
  all.grp.names <- sort(unique(all.grp.names))
  
  ## put \hphantom around the grouping variable name
  ## for all cases where a single grouping variable
  ## is used multiple times  (ie random intercept + random slope)
  for (g in 1:length(all.grp.names)) {
    
    #substr(all.ranef.names[1], 1, gregexpr(pattern =' - ', all.ranef.names[1])[[1]][1] - 1)
    raneff.gr.names <- sapply(all.ranef.names, function(x) substr(x, 1, gregexpr(pattern =' - ', x)[[1]][1] - 1))
    names(raneff.gr.names) <- NULL
    
    
    name.idx <- which(!is.na(match(raneff.gr.names, all.grp.names[g])))
    if (length(name.idx) > 1) {
      name.idx <- name.idx[-1]
      for (gg in 1:length(name.idx)) {
        all.ranef.names[name.idx[gg]] <- 
          gsub(all.grp.names[g], paste0("\\\\hphantom{", all.grp.names[g], "}"), 
               all.ranef.names[name.idx[gg]])
      }
    }
  }
  
  for (mm in 1:length(all.ranef.names)) {
    
    reps.models <- std.devs.models <- rep(" - ", n.models)
    name.raneff <- all.ranef.names[mm]
    for (m in 1:n.models) {
      
      # return the VarCorr matrix for the current model
      vc <- VarCorr(models[[m]])
      
      # get the vector of the standard deviations for
      # all random effects and the number of random 
      # effect terms for each rand eff variable
      reStdDev <- c(lapply(vc, attr, "stddev"))
      repLens <- lengths(reStdDev)
      all.std.dev <- unlist(reStdDev)
      
      # number of replications for each rand eff
      num.reps <- sapply(ranef(models[[m]]), nrow)
      
      # names of each rand effect group variable
      grp.names <- names(vc)
      cur.names <- rep(NA, sum(repLens))
      name.names <- unlist(lapply(vc, colnames))
      
      k <- 0
      for (g in 1:length(grp.names)) {
        for (gs in 1:repLens[g]) {
          k <- k + 1
          
          # get the name of each random effect
          # cur.name = Group Var - (Intercept) for the first and
          # cur.name =           - (Rand Slope) for subsequent
          # rand effects for the same grouping var
          if (gs == 1) {
            cur.names[k] <- paste(grp.names[g], name.names[k], sep = " - ")
          } else {
            
            # we need #hphantom in order to make the alignment 
            # of the random effect names work out nicely
            cur.names[k] <- paste0(paste0("\\hphantom{", grp.names[g], "}"), 
                                   " - ", name.names[k])
            
          }
        }
      }
      
      # get indices of the current random effects
      # with respect to all possible random effects
      # from all models
      ranef.idx <- match(all.ranef.names, cur.names)
      
      ## the next two lines are unnecessary. 
      reStdDev <- c(lapply(vc, attr, "stddev"))
      repLens <- lengths(reStdDev)[ranef.idx[mm]]
      
      # if the current random effect is in the current model
      # then store its name and standard deviation
      if (!is.na(ranef.idx[mm])) {
        reps.models[m] <- num.reps[ranef.idx[mm]]
        std.devs.models[m] <- prettyNum(all.std.dev[ranef.idx[mm]], big.mark=digit.separator, 
                                        scientific=FALSE, format = "f", nsmall = digits, digits = digits)
      }
    }
    
    
    number.reps <- paste0("\\# of ", name.raneff, " & ", paste(reps.models, collapse = " & "), "\\\\")
    
    stddevs <- paste0(name.raneff, " & ", paste(std.devs.models, collapse = " & "), "\\\\")
    
    #added.row.counter <- added.row.counter + 1L
    #Tables <- insertrow(Tables, number.reps, r + added.row.counter)
    
    added.row.counter <- added.row.counter + 1L
    Tables <- insertrow(Tables, stddevs,     r + added.row.counter)
    
    #if (mm < length(all.ranef.names)) {
    added.row.counter <- added.row.counter + 1L
    Tables <- insertrow(Tables, newline, r + added.row.counter)
    #}
  }
  
  std.devs.models <- rep("", n.models)
  name.raneff <- "Residual"
  for (m in 1:n.models) {
    vc <- VarCorr(models[[m]])
    
    std.devs.models[m] <- prettyNum(attr(vc, "sc"), big.mark=digit.separator, 
                                    scientific=FALSE, format = "f", nsmall = digits, digits = digits)
    
    
    
    
  }
  stddevs <- paste0(name.raneff, " Standard Deviation & ", paste(std.devs.models, collapse = " & "), "\\\\")
  added.row.counter <- added.row.counter + 1L
  Tables <- insertrow(Tables, stddevs,     r + added.row.counter)
  
  cat(paste(Tables$Tables, collapse = "\n"))
}



insertrow <- function(df, new.row, row.idx) {
  df[seq(row.idx+1, nrow(df)+1),] <- df[seq(row.idx, nrow(df)),]
  df[row.idx, ] <- new.row
  df
}