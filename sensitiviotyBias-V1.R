#Data treatments for detecting a bias from sensitivity to feedbacks
# Author: Guillaume Deffuant

library(data.table)
library(plot.matrix)

# creading data.table from file "dataExpe.csv" to data.table "dataExpe"
# dataExpe = fread("dataExpe.csv")
# downloading the data.table of results for various sets 
# allResults = fread("allResults.csv")

#-----------------------------
#computes the regressions providing approximate sensitivties for set dtli (data table)
# slope c, intercept d, sensitivity to all fbs
# slope cp, intercept dp, sensitivity to positive fbs
# slope cn, intercept dn, sensitivity to negative fbs
# enh : self-enhancement bias
# tot : total bias
# sb : sensitivity bias
# sbt : theoretical sensitivity bias
# returns all these values in a list
# threshold is a limit size of the set below which returned values are NA
valuesForSet = function(dtl, threshold){
  
  dtli = dtl[numFeedback > 0]
  if (missing(threshold)) threshold = 50
  if(length(dtli[,id]) < threshold){
    c = NA
    d = NA
    pvc = NA
    pvd = NA
  }
  else{
    lmod = summary(lm(absChgmtEvalRel ~atm1Div100, data = dtli))
    c = coef(lmod)["atm1Div100", "Estimate"]
    d = coef(lmod)["(Intercept)", "Estimate"] 
    pvc = coef(lmod)["atm1Div100", "Pr(>|t|)"]
    pvd = coef(lmod)["(Intercept)", "Pr(>|t|)"]
  }
  if(length(dtli[chgmtFeedback > 0,id]) < threshold){
    cp = NA
    dp = NA
    pvcp = NA
    pvdp = NA
  }
  else{
    lmod = summary(lm(absChgmtEvalRel ~atm1Div100, data = dtli[chgmtFeedback > 0]))
    cp = coef(lmod)["atm1Div100", "Estimate"]
    dp = coef(lmod)["(Intercept)", "Estimate"] 
    pvcp = coef(lmod)["atm1Div100", "Pr(>|t|)"]
    pvdp = coef(lmod)["(Intercept)", "Pr(>|t|)"]
  }
  if(length(dtli[chgmtFeedback < 0,id]) < threshold){
    cn = NA
    dn = NA
    pvcn = NA
    pvdn = NA
  }
  else{
    lmod = summary(lm(absChgmtEvalRel ~atm1Div100, data = dtli[chgmtFeedback < 0]))
    cn = coef(lmod)["atm1Div100", "Estimate"]
    dn = coef(lmod)["(Intercept)", "Estimate"] 
    pvcn = coef(lmod)["atm1Div100", "Pr(>|t|)"]
    pvdn = coef(lmod)["(Intercept)", "Pr(>|t|)"]
  }
  enh = mean(cp*dtli[,atm1Div100] + dp - cn*dtli[,atm1Div100] - dn)*100
  enh = round(enh, digits = 2)
  tot = mean(dtli[,a4ma0]/dtli[,abs(chgmtFeedback)])*50
  tot = round(tot, digits = 2)
  cm = (cp + cn) / 2
  dm = (dp + dn) / 2
  sbt = mean(dtli[,- cm * (cm*atm1Div100+dm)*abs(chgmtFeedback)])
  sbt2 = mean(dtli[,- cm * (cm*atm1Div100+dm)*abs(chgmtFeedback^3)])/100
  return(list('c' = c, 'd' = d, 'pvc' = pvc, 'pvd'= pvd, 'cp'= cp, 'dp' = dp, 'pvcp' = pvcp, 'pvdp' = pvdp, 'cn'= cn, 'dn' = dn, 'pvcn' = pvcn, 'pvdn' = pvdn, 'enh' = enh, 'tot' = tot, 'sbt' = sbt, 'sbt2' = sbt2))
}

#--------------------------------
# visualises the set and the sensitivities
# pres in (1, 2, 3, 4) determines the measures that are displayed
# example:
# visuForSet(dataExpe[croyance_groupe >= 7], pres = 3)
visuForSet = function(dte, threshold, pres){
  
  if (missing(threshold)) threshold = 50
  if (missing(pres)) pres = 2

  dtl = dte[numFeedback > 0 & feedbacksEquilibres > 0 & age > 15]
  
  lv = valuesForSet(dtl, threshold)
  
  atMin = min(dtl[,atm1Div100])
  atMax = max(dtl[,atm1Div100])
  par(mar=c(2.1, 2.1, 0.2, 0.1), cex = 1.5)
  plot(10, xlim = c(atMin, atMax), ylim = c(0,1), xlab = "Self-evaluation", ylab ="Self-evaluation change")
  
  yp1 = lv$cp + lv$dp
  yn1 = lv$cn + lv$dn
  if (pres >= 3){
    if ((lv$dp - lv$dn)*(lv$cp-lv$cn+lv$dp-lv$dn) < 0){
      xi = (lv$dn - lv$dp) / (lv$cp - lv$cn)
      yi = lv$cp * xi + lv$dp
      if (lv$dp > lv$dn){
        col1 = "orange"
        col2 = "lightblue"
      } 
      else{
        col1 = "lightblue"
        col2 = "orange"
      }
      polygon(c(0, xi, 0, 0), c(lv$dp, yi, lv$dn, lv$dp), col = col1, border = NULL)
      polygon(c(1, xi, 1, 1), c(yp1, yi, yn1, yp1), col = col2, border = NULL)
    }
    else{
      if (lv$dp > lv$dn) col1 = "orange"
      else col1 = "lightblue"
      polygon(c(0, 0, 1, 1, 0), c(lv$dp, lv$dn, yn1, yp1, lv$dp), col = col1, border = "black")
    }
  }
  if (pres >= 2){
    lines(c(0, 1), c(lv$dp, yp1), col = "red", lwd = 3)
    lines(c(0, 1), c(lv$dn, yn1), col = "blue", lwd = 3)
  }
  points(dtl[chgmtFeedback > 0, atm1Div100], dtl[chgmtFeedback > 0, absChgmtEvalRel], col = "red", lwd = 1)
  points(dtl[chgmtFeedback < 0, atm1Div100], dtl[chgmtFeedback < 0, absChgmtEvalRel], col = "blue", lwd = 1, pch = 3)
  if(pres >= 2){
    legend(atMin+0.1, 1, c(paste("N:", length(dtl[,id]), "  c:", round(lv$c, digits = 2), codeSignf(lv$pvc) ),paste("cp:", round(lv$cp, digits = 2), codeSignf(lv$pvcp), "  cn:", round(lv$cn, digits = 2), codeSignf(lv$pvcn))), box.col = "black", adj = 0.05)
    #print(paste("incr: ", incr))
  }
  if (pres == 3){
    legend(atMin+0.2, 0.22, c(paste("E: ", round(lv$enh, digits = 2), "  S:", round(lv$tot - lv$enh, digits = 2)), paste("T:", round(lv$tot, digits = 2), "  S':", round(lv$sbt, digits = 2))), box.col = "black", adj = 0.05)
    #  legend(atMin+0.2, 0.22, paste("E: ", se, "  S:", dt), box.col = "black", adj = 0.05)
  }
  
  if (pres == 4){
    enh = round(lv$enh, digits = 2)
    sbt = round(lv$sbt, digits = 2)
    legend(atMin+0.2, 0.22, c(paste("E: ", enh), paste("S':", sbt)), box.col = "black", adj = 0.05)
    #  legend(atMin+0.2, 0.22, paste("E: ", se, "  S:", dt), box.col = "black", adj = 0.05)
  }
  
  lines(c(0, 1), c(lv$d, lv$c+lv$d), col = "black", lwd = 3)
  
  if(pres == 1){
    legend(atMin+0.15, 1, paste("N:", length(dtl[,id]),"  c:", round(lv$c, digits = 2), codeSignf(lv$pvc)), box.col = "black", adj = 0.05)
  }
}


#-------------------
# code for significance from p value for latex
codeSignfTex = function(pv){
  
  if(is.na(pv)) return("")
  if(pv < 0.001) return("^{***}")
  if(pv < 0.01) return("^{**}")
  if(pv < 0.05) return("^*")
  if(pv < 0.1) return("\\hspace{0.1 cm}.")
  #if(pv < 0.1) return(" .")
  return("")
}

#-------------------
# code for significance from p value display in R plot
codeSignf = function(pv){
  
  if(is.na(pv)) return("")
  if(pv < 0.001) return("***")
  if(pv < 0.01) return("**")
  if(pv < 0.05) return("*")
  if(pv < 0.1) return(" .")
  return("")
}



#---------------------------
#bootsrap on the measures of bias
# Cas where the 4 time steps are present
# bootstrap on individuals
bootIndValuesForSet = function(dt, nRep, threshold){
  
  nb = length(dt[,id])/4
  lenh = NULL
  ltot = NULL
  lsbt = NULL
  lsb = NULL
  ldiff = NULL
  for (i in (1:nRep)){
    indsId = sample(nb, replace = T)
    inds = NULL
    for (j in (1:nb)) inds = c(inds, (indsId[j] * 4) - (0:3))
    lv = measuresForSet(dt[inds,], threshold)
    lenh = c(lenh, lv$enh)
    ltot = c(ltot, lv$tot)
    lsbt = c(lsbt, lv$sbt)
    lsb = c(lsb, lv$sb)
    ldiff = c(ldiff, lv$sb - lv$sbt)
  }
  return(list('enhm' = mean(lenh), 'enhsd'= sd(lenh), 'totm' = mean(ltot), 'totsd' = sd(ltot), 'sbtm' = mean(lsbt), 'sbtsd' = sd(lsbt), 'sbm' = mean(lsb), 'sbsd' = sd(lsb), 'diffm' = mean(ldiff), 'diffsd' = sd(ldiff)))
}

#---------------------------
#bootsrap on the measures of bias
# case where not all time steps are present
bootValuesForSet = function(dt, nRep, threshold){
  
  nb = length(dt[,id])
  lenh = NULL
  ltot = NULL
  lsbt = NULL
  lsb = NULL
  ldiff = NULL
  for (i in (1:nRep)){
    inds = sample(nb, replace = T)
    lv = measuresForSet(dt[inds,], threshold)
    lenh = c(lenh, lv$enh)
    ltot = c(ltot, lv$tot)
    lsbt = c(lsbt, lv$sbt)
    lsb = c(lsb, lv$sb)
    ldiff = c(ldiff, lv$sb - lv$sbt)
  }
  return(list('enhm' = mean(lenh), 'enhsd'= sd(lenh), 'totm' = mean(ltot), 'totsd' = sd(ltot), 'sbtm' = mean(lsbt), 'sbtsd' = sd(lsbt), 'sbm' = mean(lsb), 'sbsd' = sd(lsb), 'diffm' = mean(ldiff), 'diffsd' = sd(ldiff)))
}


#--------------------------------
# Returns a data.table with the sensitivities and all measures and N
# for 7938 sets (see in the beginning the values of interviewTime, trust, anchor,...)
# the bootstrap is activated if nRep > 0
# File "allResults.csv" is the result of this function for nRep = 200
tabAllMeasures = function(dte, threshold, ageMin, nRep){
  
  if(missing(threshold)) threshold = 0 
  if (missing(ageMin)) ageMin = 15
  if (missing(nRep)) nRep = 0
  tints = list(c(1,1), c(2,2), c(3,3), c(4,4), c(1,2), c(1,3), c(1,4))
  #intTimeMins = c(0, 240)
  #trustInts = list(c(0,10), c(0,3), c(0,4), c(0,5), c(6,10), c(7,10), c(8,10), c(9,10)) 
  trustInts = list(c(0,10), c(0,5), c(0,6), c(6,10), c(7,10), c(8,10), c(9,10)) 
  anchorInts = list(c(0,100), c(0,50), c(50,100))
  seInts = list(c(0,5), c(0,3), c(3.01,5))
  genderVals = c("all", "Un homme", "Une femme")
  scaleVals = c("all", 0, 1)
  #nst = length(tints)*length(intTimeMins)*length(trustInts)*length(anchorInts)*length(seInts)*length(genderVals)*length(scaleVals)
  nst = length(tints)*length(trustInts)*length(anchorInts)*length(seInts)*length(genderVals)*length(scaleVals)
  trustCol = NULL
  anchorCol = NULL
  seCol = NULL
  genderCol = NULL
  scaleCol = NULL
  NCol = NULL
  cCol = NULL
  pvcCol = NULL
  dCol = NULL
  pvdCol = NULL
  cpCol = NULL
  pvcpCol = NULL
  dpCol = NULL
  pvdpCol = NULL
  cnCol = NULL
  pvcnCol = NULL
  dnCol = NULL
  pvdnCol = NULL
  enhCol = NULL
  totCol = NULL
  sbtCol = NULL
  sbt2Col = NULL
  tCol = NULL
  intTimeCol = NULL
  if (nRep > 0){
    enhmCol = NULL
    enhsdCol = NULL
    totmCol = NULL
    totsdCol = NULL
    sbmCol = NULL
    sbsdCol = NULL
    sbtmCol = NULL
    sbtsdCol = NULL
    diffmCol = NULL
    diffsdCol = NULL
  }
  dtl = dte[feedbacksEquilibres > 0 & age >ageMin & numFeedback > 0]
  currs = 0
  for (itrust in(1:length(trustInts))){
    for (ianchor in(1:length(anchorInts)) ){
      for (ise in (1:length(seInts))){
        for (igender in (1:length(genderVals))){
          for(iscale in (1:length(scaleVals))){
            for (it in (1:length(tints))){
              # for(intTime in (1:length(intTimeMins))){
              currs = currs + 1
              cat("\r ", currs, " over  ", nst)
              dtli = dtl[croyance_groupe >= trustInts[[itrust]][1] & croyance_groupe <= trustInts[[itrust]][2]]
              dtli = dtli[ancrage >= anchorInts[[ianchor]][1] & ancrage <= anchorInts[[ianchor]][2]]
              dtli = dtli[selfEsteem >= seInts[[ise]][1] & selfEsteem <= seInts[[ise]][2]]
              if(igender > 1) dtli = dtli[genre == genderVals[igender]]
              if(iscale > 1) dtli = dtli[echelleCroissante == scaleVals[iscale]]
              dtli = dtli[numFeedback >= tints[[it]][1] & numFeedback <=  tints[[it]][2]]
              # dtli = dtli[interviewtime > intTimeMins[intTime]]
              trustCol = c(trustCol, paste0("[",trustInts[[itrust]][1],", ", trustInts[[itrust]][2], "]"))
              anchorCol =  c(anchorCol, paste0("[",anchorInts[[ianchor]][1],", ", anchorInts[[ianchor]][2], "]"))
              seCol = c(seCol, paste0("[",seInts[[ise]][1],", ", seInts[[ise]][2], "]"))
              genderCol = c(genderCol, genderVals[igender])
              scaleCol = c(scaleCol, scaleVals[iscale])
              NCol = c(NCol, length(dtli[,id]))
              tCol = c(tCol, paste0(tints[[it]][1],":",tints[[it]][2]))
              #intTimeCol = c(intTimeCol, intTimeMins[intTime])
              lv = valuesForSet(dtli, threshold)
              cCol = c(cCol, lv$c)
              dCol = c(dCol, lv$d)
              pvcCol = c(pvcCol, lv$pvc)
              cpCol = c(cpCol, lv$cp)
              dpCol = c(dpCol, lv$dp)
              pvcpCol = c(pvcpCol, lv$pvcp)
              pvdpCol = c(pvdpCol, lv$pvdp)
              cnCol = c(cnCol, lv$cn)
              dnCol = c(dnCol, lv$dn)
              pvcnCol = c(pvcnCol, lv$pvcn)
              pvdnCol = c(pvdnCol, lv$pvdn)
              totCol = c(totCol, lv$tot)
              sbtCol = c(sbtCol, lv$sbt)
              sbt2Col = c(sbt2Col, lv$sbt2)
              if (nRep > 0){
                if (tints[[it]][1] == 1 & tints[[it]][2] == 4) lv = bootIndValuesForSet(dtli, nRep, threshold)
                else lv = bootValuesForSet(dtli, nRep, threshold)
                enhmCol = c(enhmCol, lv$enhm)
                enhsdCol = c(enhsdCol, lv$enhsd)
                totmCol = c(totmCol, lv$totm)
                totsdCol = c(totsdCol, lv$totsd)
                sbtmCol = c(sbtmCol, lv$sbtm)
                sbtsdCol = c(sbtsdCol, lv$sbtsd)
                sbmCol = c(sbmCol, lv$sbm)
                sbsdCol = c(sbsdCol, lv$sbsd)
                diffmCol = c(diffmCol, lv$diffm)
                diffsdCol = c(diffsdCol, lv$diffsd)
              }
            }
          }
        }
      }
    }
  }
  # }
  if (nRep == 0) res = data.table(
    trust = trustCol,
    anchor = anchorCol,
    se = seCol,
    gender = genderCol,
    scale = scaleCol,
    t = tCol,
    # intTime = intTimeCol,
    N = NCol,
    c = cCol,
    pvc = pvcCol,
    d = dCol,
    pvd = pvdCol,
    cp = cpCol,
    pvcp = pvcpCol,
    dp = dpCol,
    pvdp = pvdpCol,
    cn = cnCol,
    pvcn = pvcnCol,
    dn = dnCol,
    pvdn = pvdnCol,
    enh = enhCol,
    tot = totCol,
    sbt = sbtCol,
    sbt2 = sbt2Col
  )
  else res = data.table(
    trust = trustCol,
    anchor = anchorCol,
    se = seCol,
    gender = genderCol,
    scale = scaleCol,
    t = tCol,
    # intTime = intTimeCol,
    N = NCol,
    c = cCol,
    pvc = pvcCol,
    d = dCol,
    pvd = pvdCol,
    cp = cpCol,
    pvcp = pvcpCol,
    dp = dpCol,
    pvdp = pvdpCol,
    cn = cnCol,
    pvcn = pvcnCol,
    dn = dnCol,
    pvdn = pvdnCol,
    enh = enhCol,
    tot = totCol,
    sbt = sbtCol,
    sbt2 = sbt2Col,
    enhm = enhmCol,
    enhsd = enhsdCol,
    totm = totmCol,
    totsd = totsdCol,
    sbm = sbmCol,
    sbsd = sbsdCol,
    sbtm = sbtmCol,
    sbtsd = sbtsdCol,
    diffm = diffmCol,
    diffsd = diffsdCol
  )
}
#----------------------------------
#Hierarchical linear models

#--------------------------------
# hlm1 level1
# returns the slope and intercept of sensitivity for each intididual
# with significance
# and other variables (anchor, meana, age...) for the second level regressions
hlm1level1 = function(dte){
  
  dtIndiv = dte[numFeedback == 0]
  dtl = dte[numFeedback > 0]
  ids = dtIndiv[,id]
  dtl[,atm1 := atm1/100]
  coefs = NULL
  ints = NULL
  pvcs = NULL
  pvis = NULL
  for (ind in ids){
    lmodp = summary(lm(absChgmtEvalRel ~atm1Div100, data = dtl[id == ind]))
    if (length(coef(lmodp)[,"Estimate"]) < 2)  {
        coefs = c(coefs, NA)
        ints = c(ints, NA)
        pvcs = c(pvcs, NA)
        pvis = c(pvis, NA)
    }
    else{
        coefs = c(coefs, coef(lmodp)["atm1Div100", "Estimate"])
        ints = c(ints, coef(lmodp)["(Intercept)", "Estimate"])
        pvcs = c(pvcs, coef(lmodp)["atm1Div100", "Pr(>|t|)"])
        pvis = c(pvis, coef(lmodp)["(Intercept)", "Pr(>|t|)"])
    }
  }
 age = dtIndiv[, age]
 meana = dtIndiv[,meana]
 anchor = dtIndiv[,ancrage]
 se = dtIndiv[,selfEsteem]
 dat = list("coefs" = coefs, "ints" = ints, "pvcs" = pvcs, "pvis" = pvis, "age" = age, "meana" = meana, "anchor"= anchor, "se" = se)
}


#Example of second level regression on l1data returned by hlm1level1:
# l1data = hlm1leve1(dataExpe)
# summary(lm(coefs ~ meana, data = l1data))



#--------------------------------
# hierachical model 2 (looking for pattern in time)
# Returns the slope and intercept of second level regressions 
# slope ~ time or
# intercept ~ time
# with time = 1:4
hlm2 = function(dte){
  
  threshold = 50
  dtl = dte[feedbacksEquilibres > 0 & numFeedback > 0]
  ints = rep(0, 4)
  coefs = rep(0, 4)
  intsp = rep(0, 4)
  coefsp = rep(0, 4)
  intsn = rep(0, 4)
  coefsn = rep(0, 4)
  for (tt in 1:4){
    lv = valuesForSet(dtl[numFeedback == tt], 50)
    ints[tt] = lv$d
    coefs[tt] = lv$c
    intsp[tt] = lv$dp
    coefsp[tt] = lv$cp
    intsn[tt] = lv$dn
    coefsn[tt] = lv$cn
  }
  regs = list("time" = 1:4, "ints" = ints, "coefs" = coefs, "intsp" = intsp, "coefsp" = coefsp, "intsn" = intsn, "coefsn" = coefsn)
  #  print(regs)
  lms = list(0, 0, 0, 0, 0, 0)
  names(lms) = c("ints", "coefs","intsp", "coefsp", "intsn", "coefsn")
  for(tt in 1:6){
    lms[[tt]] = summary(lm(regs[[tt+1]] ~time, data = regs))
    print(paste0(names(lms)[tt],":  ",round(coef(lms[[tt]])["time", "Estimate"], digit = 2), codeSignf(coef(lms[[tt]])["time", "Pr(>|t|)"])))
  }
}



