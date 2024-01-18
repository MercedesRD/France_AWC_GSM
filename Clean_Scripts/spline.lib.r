# Name: spline.lib.r
# Description: tools for using the spline 
# O/S: any
# Date: 03/07/2013
# Author: Marine Lacoste
# Company: INRA

# N.B. : Fonctions écrites par Manuel Martin, Rossanno Ciampalini et Marine Lacoste




## Liste des fonctions du répertoire
###########################################################################################################################################################

## Fonctions pour le calcul des splines sur un profil
#...................................................
# splineCostFunc  : fonction de coût qui ajuste un spline sur un profil pour un lambda donné et renvoie le tmse
# getOptimumLambda : uses the R optimize function to get the ldba parameter (firest parameter) yielding the minimum value
# fitSpline : fit a spline fonction to a soil profile data
# predSplineSingle : fonction permettant de prédire les valeurs d'une propriété de sol pour une profondeur donnée
# predSpline : fonction permettant de prédire les valeurs d'une propriété de sol pour un profil et pour une série de profondeurs déterminées
#               (fonction basée sur predSplineSingle)
# intgSpline : predictions for the defined horizons from a depth fonction (integration of a spline function)


## Fonctions pour le calcul des splines sur plusieurs profils
#............................................................
# pred.spl.multiP1 : fonction permettant de prédire les valeurs d'une propriété de sol continue pour des horizons standards et pour de multiples profils.
# pred.spl.multiP2 : fonction permettant de prédire les valeurs d'une propriété de sol continue pour des horizons standards et pour de multiples profils.
# pred.spl.multiPV : fonction permettant de prédire les valeurs d'une propriété de sol continue pour des horizons standards pour plusieurs variables et pour de multiples profils.

# Différece entre pred.spl.multiP1 et pred.spl.multiP2 : structure des données en sortie

## Fonction pour la préparation des horizons de surface pour les profils labourés
#................................................................................
# prep.hz : préparation des horizons pour spline : ajout de deux horizons d'épaisseur xd aux extrémités des horizons labourés (fonction pour un profil)
# prep.hz.multi : préparation des horizons pour spline : ajout de deux horizons d'épaisseur xd aux extrémités des horizons labourés (fonction pour plusieurs profils)

## Fonctions pour la représentation graphique des résultats
#..........................................................
# max.min.var : function to determine the max and min of the analyte for plotting
# plot_raw : fonction permettant de représenter graphiquement les valeurs observées sur un profil 
# plot_raw.sp : fonction permettant de représenter graphiquement les valeurs observées sur un profil, plus les valeurs prédites
#               par le spline  (par centimètres uniquement)
# plotAll : fonction permettant de représenter graphiquement les valeurs observées sur un profil, plus les valeurs prédites
#               par le spline  (par centimètres et pour des horizons standards)
# plotSeq : fonction permettant de représenter graphiquement les résultats des splines pour plusieurs profils 
#               et d'exporter le résultat sous forme d'un pdf

## Autres fonctions
#..........................................................
# conc.champ : fonction permettant rajouter une colonne à une dataframe, contenat
# la concatenation des champs provenant de plusieurs colonnes
###########################################################################################################################################################

#' Prediction unique sur un spline
#'
#' Fonction permettant de prédire les valeurs d'une propriété de sol pour un profil et pour une série de profondeurs déterminées
#' @param z : une profondeur unique
#' @param spl : modèle Spline de type liste dont les paramètres sont f, b0, b1, delta, d (sortie de la fonction fitspline)
#' @return un vecteur de longueur 1 donnant la valeur prédite pour une profondeur donnée
#' @export 
#' @keywords spline
predSplineSingle <- function(z, spl){
  
  ## tests which intervals of the spline z belongs to
  z <- min(z, max(spl$lower) - 1e-6)
  ## rescale the depth
  u <- spl$upper
  v <- spl$lower
  
  z <- min(z, max(v) - 1e-6)
  b1 <- spl$b1
  b0 <- spl$b0
  gamm <- spl$gamm
  alfa <- spl$alfa
  
  hz <- which(z < spl$lower & z >= spl$upper)
  if (length(hz) > 0){#belongs to an observed horizon
    res <- alfa[hz]+b0[hz]*(z-u[hz])+gamm[hz]*(z-u[hz])^2
  } else { # lies between two Observed horizons
    hz <- which(z >= c(0, spl$lower) & z < c(spl$upper, 0)) - 1
    phhz <- alfa[hz+1] - b1[hz] * (u[hz+1] - v[hz])
    res <- phhz + b1[hz] %*% (z - v[hz])
  }
  
  ## performs the prediction
  return( as.numeric(res) )
}

#' Predictions multiples sur un spline
#' 
#' Fonction permettant de prédire les valeurs d'une propriété de sol pour un profil et pour une série de profondeurs déterminées
#' @param zs : un vecteur de profondeurs
#' @param spl : modèle Spline de type liste dont les paramètres sont f, b0, b1, delta, d (sortie de la fonction fitspline)
#' @return un vecteur de longueur de zs donnant les valeur prédite pour les
#' profondeurs données
#' @export
#' @keywords spline
predSpline <- function(zs, spl){
  ### DESCRIPTION : fonction permettant de prédire les valeurs d'une propriété de sol pour un profil et pour une série de profondeurs déterminées
  ### INPUTS :
  # - zs : un vecteur de profondeurs
  # - spl : modèle Spline de type liste dont les paramètres sont f, b0, b1, delta, d (sortie de la fonction fitspline)
  ### OUTPUTS :
  # - un vecteur de la longueur de zs donnat les valeurs prédites pour les profondeurs données
  
  sapply(zs, function(x){predSplineSingle(x, spl)}, simplify = T)
}


#' fits a spline fonction to a soil profile data
#' 
#' The profil must have at least 2 horizons but these horizons need not to be
#' continuous. This make the procedure applicable for a dataset describing a set of
#' soil samples, for wich the upper and lower boundaries are known, and with soil
#' property measurement we want to spline on. 
#' N.B. The initial scipt for this function comes frm Budiman Minasny (Sydney University)
#' @param y : soil horizon values to fit (vector, mandatory)
#' @param upper.boundary : upper limits of the described soil horizons (vector, mandatory)
#' @param lower.boundary : upper limits of the described soil horizons (vector, mandatory)
#' @param lambda : the smoothing parameter (a real, mandatory)
#' @return a list of 4 elements
#' - f : the fit values at the soil horizon boundaries (one column matrix, with one line per horizon boundary)
#' - b0 : spline parameters (one column matrix, with one line per horizon)
#' - b1 : spline parameters (one column matrix, with one line per horizon)
#' - tmse : estimate of the true mean squared error (a real)
#' @keywords profile spline depth soil property interpolation
#' @examples
#' test <- data.frame(width = c(1, 12, 5, 4), clay = c(2, 30, 50, 20), upper = c(0, 1,
#' 			13, 18), lower = c(1, 13, 18, 22))
#' splin <- fitSpline(test$clay, test$upper, test$lower, 0.0001)
#' plot(test$clay, -(test$upper + test$lower) / 2, ylim = -range(test$upper,
#' test$lower), xlim = c(0, max(test$clay) * 1.5))
#' lines(predSpline(0:max(test$lower), splin), -(0:max(test$lower)))
#' @export
#' @references Malone BP, McBratney AB, Minasny B, Laslett GM: Mapping continuous depth functions of soil carbon storage and available water capacity. Geoderma 2009, 154:138-152.
fitSpline <-  function(y,upper.boundary,lower.boundary,lambda){
 
  s <- 0.05*sd(y)       # 5% of the standard deviation of the target attribute
  u <- upper.boundary
  v <- lower.boundary
  x <- c(0:max(v))
  
  # n is the number of observations
  # so n+1 is the number of interval boundaries
  n<- length(y)
  np1<- n+1
  nm1<- n-1
  nm2<- n-2
  nm3<- n-3
  
  # delta is v-u
  delta<- v-u
  # del is (u1-v0,u2-v1, ...)
  #del<- (u[2:n],u[n])-v
  
  del<-u
  del[1:n-1]<-u[2:n]
  del<-del-v
  
  #del<-u[2:n,n]-v
  
  # create the (n-1)x(n-1) matrix r
  # first create r with 1's on the diagonal and upper diagonal, and
  # 0's elsewhere
  vr<- c(1,1,rep(0,nm2))
  vr<- c(rep(vr,nm2),1)
  r<- matrix(vr,nm1,nm1,byrow=T)
  
  # then create a diagonal matrix d2 of differences to premultiply
  # the current r
  d2<- matrix(0,nm1,nm1)
  diag(d2)<- delta[2:n]
  # then premultiply and add the transpose; this gives half of r
  r<- d2%*%r
  r<- r+t(r)
  # then create a new diagonal matrix fo differences to add to the diagonal
  d1<- matrix(0,nm1,nm1)
  diag(d1)<- delta[1:nm1]
  
  d3<- matrix(0,nm1,nm1)
  diag(d3)<- del[1:nm1]
  
  r<- r+2*d1+6*d3       ########Different to old one
  
  # create the (n-1)xn matrix qt (transpose of q)
  vr<- c(-1,1,rep(0,nm1))
  vr<- c(rep(vr,nm2),-1,1)
  qt<- matrix(vr,nm1,n,byrow=T)
  
  
  # inverse of r
  rinv<- solve(r)
  
  # identity matrix i
  i<- diag(rep(1,n))
  
  # create the matrix coefficient z of f
  rq <- rinv %*% qt
  z <- 6 * n * lambda * t(rq) %*% qt + i
  
  # solve for the fitted horizon means
  sbar <- solve(z,y)             ############different object created
  
  # calculate the fitted values at the knots
  b<- 6*rq%*%sbar
  b0<- c(0,b)
  b1<- c(b,0)
  
  gamm<-(b1-b0)/t(2*delta);
  #Warning message:
  #  In (b1 - b0)/t(2 * delta) :
  #  longer object length is not a multiple of shorter object length
  #browser()
  alfa<-sbar-b0*t(delta)/2-gamm*t(delta)^2/3;
  
  ssq <-  sum((t(y)-sbar)^2)
  g <-  solve(z)
  ei <-  eigen(g)
  ei <-  ei$values
  df <-  n-sum(ei)
  sig2w <-  ssq/df
  # calculate the Carter and Eagleson estimate of residual variance
  dfc <-  n-2*sum(ei)+sum(ei^2)
  sig2c <-  ssq/dfc
  # calculate the estimate of the true mean squared error
  tmse <- ssq / n - 2 * s^2 * df / n + s^2
  
  ##########################################################################################################  
  # Outputs
  res  <-  list(f = NA, b0 = b0, b1 = b1, tmse = tmse, gamm = gamm, alfa = alfa, upper = u, lower = v)
  
  return(res)
  
}
test <- data.frame(width = c(1, 12, 5, 4), clay = c(2, 30, 50, 20), upper = c(0, 1,
			13, 18), lower = c(1, 13, 18, 22))
splin <- fitSpline(test$clay, test$upper, test$lower, lambda = 0.0014488)
plot(test$clay, -(test$upper + test$lower) / 2, ylim = -range(test$upper,
test$lower), xlim = c(0, max(test$clay) * 1.5))
lines(predSpline(0:max(test$lower), splin), -(0:max(test$lower)))



#' Integrates Spline functions
#'
#' Intregates prediction for the defined horizons from a depth fonction (integration of a spline function)
#' @param std : limits of the horizon boundaries (vector, mandatory)
#' @param spl : modèle Spline de type liste dont les paramètres sont f, b0, b1, delta, d (sortie de la fonction fitspline)
#' @param averaged : does one want true integral or averaged property?
#' @return  predicted value for each horizon (vector)
#' @examples
#' test <- data.frame(width = c(1, 12, 5, 4), clay = c(2, 30, 50, 20), upper = c(0, 1,
#' 			13, 18), lower = c(1, 13, 18, 22))
#' splin <- fitSpline(test$clay, test$upper, test$lower, lambda = 0.0014488)
#' intgSpline(c(0, 30), splin)
#' @export
intgSpline <- function(std, spl, averaged = T){
  # DESCRIPTION : predictions for the defined horizons from a depth fonction (integration of a spline function)
  # INPUTS
  # - std : limits of the horizon boundaries (vector, mandatory)
  # - spl : spline model produce by the fitSpline fonction (list, mandatory)
  # OUTPUTS : 
  # - yint : predicted value for each horizon (vector)
  
  # nombre d'horizons
  nd <-  length(std)-1
  yint <- rep(NA,nd)
  
  for (i in 1:nd) {
    xd1 <-  std[i]
    xd2 <-  std[i+1] 
    ## AVERAGED? does one want true integral or average property?
    if (averaged)
	width <- xd2-xd1 else
	width <- 1
    yint[i]<- (integrate(predSpline,xd1,xd2,spl)$value)/(width)
  }
  return(yint)
} 

splineCostFunc <-  function(lbda, Y, U, V){
  ## DESCRIPTION : Fonction de Coût qui ajuste un spline sur un profil pour un lambda donné et renvoie le tmse
  ## INPUTS
  # - Paramètres : lambda (main spline's parameter)
  # - Y : Response variable values
  # - U : upper boundary of the soil layers matching the response values
  # - V : lower boundary of the soil layers matching the response values
  ## OUTPUT
  # - the TMSE value
  
  splin <- fitSpline(Y, U, V, lbda)
  
  return(splin$tmse)
}

#' Uses the R optimize funciton to get the ldba parameter (firest parameter) yielding the minimum value
#'                of the splineCostFunction, that is minimum tmse values, for a given profile.
#' @param Y : Response variable values
#' @param U : upper boundary of the soil layers matching the response values
#' @param V : lower boundary of the soil layers matching the response values
#' @param minMax : vecteur de longueur 2 donnant les valeurs min et max des lambda à tester
#' @return the value of the optimum lambda parameter
#' @keywords profile spline optimization
#' @examples 
#' test <- data.frame(width = c(0, 12, 5, 4), clay = c(2, 30, 50, 20), upper = c(0, 0,
#' 			12, 17), lower = c(0, 12, 17, 21))
#' getOptimumLambda(test$clay, test$upper, test$lower)
#' @export
getOptimumLambda <- function(Y, U, V, minMax = c(1e-10, 100)){
  ## DESCRIPTION : Uses the R optimize funciton to get the ldba parameter (firest parameter) yielding the minimum value
  #                of the splineCostFunction, that is minimum tmse values, for a given profile.
  ## INPUTS : 
  # - Y : Response variable values
  # - U : upper boundary of the soil layers matching the response values
  # - V : lower boundary of the soil layers matching the response values
  # - minMax : vecteur de longueur 2 donnant les valeurs min et max des lambda à tester
  ## OUTPUT
  #  - the value of the optimum lambda parameter
  
  outM <- optimize(f = splineCostFunc, minMax, Y = Y, U = U, V = V)
  return(outM$minimum)
}


pred.spl.multiP1 <- function(data, std.depth, champ.id, champ.prof.sup, champ.prof.inf, champ.variable){
  ## DESCRIPTION : Fonction permettant de prédire les valeurs d'une propriété de sol continue pour des horizons standards et pour de multiples profils.
  ## INPUTS :
  # - data : data.frame comportant les profils de sols
  # - champ.id : nom de la colonne contenant les numéro de profil
  # - champ.prof.sup : nom de la colonne contenant les limites supérieures de profil
  # - champ.prof.inf : nom de la colonne contenant les limites inférieures de profil
  # - champ.variable : nom de la colonne contenant les valeurs de la propriété de sol à prédire
  ## OUPUTS : Liste de 2 éléments
  # - 1er élément : data.frame contenant les valeurs estimées par les splines pour chaque horizon standard (un profil par ligne)
  #       -->      2 colonnes systématiques : Soil.ID, lambda utilisé
  #       -->      une colonne par horizon standard
  #       -->      une colonne contenant la profondeur max des horizon standard 
  # - 2nd élément : data.frame contenant les valeurs estimées par les splines pour cm sur les horizons standards (un profil par colonne)
  
  # N.B. : cette fonction est à utiliser pour représenter graphiquement le résultat des splines avec la fonction plotSeq
  
  data.spl <- subset(data, select = c(champ.id, champ.prof.sup, champ.prof.inf, champ.variable))
  liste <- unique(data.spl[,1])
  liste.spl <- vector()
  liste.nspl <- vector()
  res.hz <- matrix()
  res.cm <- matrix()
  
  # Initialisation de la barre de progression  
  pb <- txtProgressBar(min = 0, max = length(liste) , style = 3)
  
  for(i in 1:length(liste)){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    # print(paste(i, liste[i]))
    prof.i <- data.spl [data.spl[,champ.id] == liste[i],]
    prof.i <- prof.i[order(prof.i[,champ.prof.sup],decreasing =F),] 
    prof.i <- prof.i[!is.na(prof.i[,champ.variable]),]
    
    if (nrow(prof.i) > 1 & min(prof.i[,champ.prof.sup]) == 0){
      upper.i <- prof.i[,champ.prof.sup]
      lower.i <- prof.i[,champ.prof.inf]
      values.i <- prof.i[,champ.variable]
      
      # sélection de meilleur lambda
      lambda <- getOptimumLambda(values.i, upper.i, lower.i, minMax = c(1e-10, 100))
      
      # Calcul du spline 
      mod.spl <- fitSpline(values.i,upper.i,lower.i,lambda)
      
      # prediction pour les horizons standards
      pred.hz <- intgSpline(std.depth, mod.spl)
      pred.hz <- c(prof.i[i,champ.id], lambda, pred.hz,200)
      # prediction par cm sur toute la profondeur des horizons standards
      pred.cm <- predSpline(seq(0, max(std.depth),1), mod.spl)
      
      # concaténation des résultats pour les horizons standards
      if (i == 1){
        res.hz <- pred.hz
      } else {
        res.hz <- rbind(res.hz, pred.hz)
      }
      # concaténation des résultats par centimètre
      if (i == 1){
        res.cm <- pred.cm
      } else {
        res.cm <- cbind(res.cm, pred.cm)
      }
      
    } else {
      # cas où le profil n'est composé que d'un seul horizon (pas de spline possible)
      liste.nspl <- c(liste.nspl, data.spl[,1])
    }
  }
  
  close(pb)  
  
  for (j in 1:(length(std.depth)-1)){
    if (j == 1){
      nom <- paste("H",std.depth[j],std.depth[j+1], sep=".")  
    } else {
      nom <- c(nom, paste("H",std.depth[j],std.depth[j+1], sep="."))
    }
  }
  
  # Mise en forme des résultats
  res.hz <- as.data.frame(res.hz)
  colnames(res.hz) <- c(champ.id, "lambda", nom, "prof.max")
  rownames(res.hz) <- seq(1, nrow(res.hz), 1)
  
  res.cm <- as.data.frame(res.cm)
  colnames(res.cm) <- liste
  rownames(res.cm) <- seq(0, max(std.depth),1)
  
  return(list(res.hz, res.cm))
}

pred.spl.multiP2 <- function(data, std.depth, champ.id, champ.prof.sup, champ.prof.inf, champ.variable){
  ## DESCRIPTION : fonction permettant de prédire les valeurs d'une propriété de sol continue pour des horizons standards et pour de multiples profils.
  ## INPUTS:
  # - data : data.frame comportant les profils de sols
  # - champ.id : nom de la colonne contenant les numéro de profil
  # - champ.prof.sup : nom de la colonne contenant les limites supérieures de profil
  # - champ.prof.inf : nom de la colonne contenant les limites inférieures de profil
  # - champ.variable : nom de la colonne contenant les valeurs de la propriété de sol à prédire
  ## OUPUTS: liste de deux éléments 
  # - 1er élément : data.frame contenant les valeurs estimées par les splines pour chaque horizon standard 
  #             4 colonnes : Soil.ID, Upper.Boundary, Lower.Boundary, SOC
  # - 2nd élément : data.frame contenant les valeurs estimées par les splines pour chaque cm. une colonne par profil
  
  data.spl <- subset(data, select = c(champ.id, champ.prof.sup, champ.prof.inf, champ.variable))
  liste <- unique(data.spl[,champ.id])
  liste.spl <- vector()
  liste.nspl <- vector()
  res.hz <- matrix()
  res.cm <- matrix()
  for (j in 1:(length(std.depth)-1)){
    if (j == 1){
      nom <- paste("H",std.depth[j],std.depth[j+1], sep=".")  
    } else {
      nom <- c(nom, paste("H",std.depth[j],std.depth[j+1], sep="."))
    }
  }
  
  # Initialisation de la barre de progression  
  pb <- txtProgressBar(min = 0, max = length(liste) , style = 3)
  
  for(i in 1:length(liste)){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    prof.i <- data.spl[data.spl[,champ.id] == liste[i],]
    prof.i <- prof.i[order(prof.i[,champ.prof.sup],decreasing =F),] 
    prof.i <- prof.i[!is.na(prof.i[,champ.variable]),]
    
    if (nrow(prof.i) > 1 & min(prof.i[,champ.prof.sup]) == 0){
      
      liste.spl <- c(liste.spl,unique(prof.i[,champ.id]))
      
      upper.i <- prof.i[,champ.prof.sup]
      lower.i <- prof.i[,champ.prof.inf]
      values.i <- prof.i[,champ.variable]
      
      # sélection de meilleur lambda
      lambda <- getOptimumLambda(values.i, upper.i, lower.i, minMax = c(1e-10, 100))
      
      # Calcul du spline  
      mod.spl <- fitSpline(values.i,upper.i,lower.i,lambda)
      
      # prediction pour les horizons standards
      pred.hz <- intgSpline(std.depth, mod.spl)
      pred.hz <- cbind(rep(prof.i[1,champ.id],length(std.depth)-1), std.depth[1:length(std.depth)-1], std.depth[2:length(std.depth)], pred.hz)
      # prediction par cm sur toute la profondeur des horizons standards
      pred.cm <- predSpline(seq(0, max(std.depth),1), mod.spl)
      
      # concaténation des résultats pour les horizons standards
      if (i == 1){
        res.hz <- pred.hz
      } else {
        res.hz <- rbind(res.hz, pred.hz)
      }
      # concaténation des résultats par centimètre
      if (i == 1){
        res.cm <- pred.cm
      } else {
        res.cm <- cbind(res.cm, pred.cm)
      }
      
    } else {
      # cas où le profil n'est composé que d'un seul horizon (pas de spline possible)
      liste.nspl <- c(liste.nspl, data.spl[,champ.id])
    }
  }
  
  close(pb)
  
  # Mise en forme des résultats
  res.hz <- as.data.frame(res.hz)
  colnames(res.hz) <- c(champ.id, champ.prof.sup, champ.prof.inf, champ.variable)
  rownames(res.hz) <- seq(1, nrow(res.hz), 1)
  
  res.cm <- as.data.frame(res.cm)
  colnames(res.cm) <- liste
  rownames(res.cm) <- seq(0, max(std.depth),1)
  
  
  return(list(res.hz, res.cm))
}

pred.spl.multiPV <- function(data, std.depth, champ.id, champ.prof.sup, champ.prof.inf, champ.variables){
  ## DESCRIPTION : fonction permettant de prédire les valeurs d'une propriété de sol continue pour des horizons standards
  ##                pour plusieurs variables et pour de multiples profils.
  ## INPUTS 
  ## - data : data.frame comportant les profils de sols
  ## - champ.id : nom de la colonne contenant les numéro de profil
  ## - champ.prof.sup : nom de la colonne contenant les limites supérieures de profil
  ## - champ.prof.inf : nom de la colonne contenant les limites inférieures de profil
  ## - champ.variables : liste des colonnes contenant les valeurs des propriétés de sol à prédire
  ## OUTPUTS 
  # - data.frame contenant les valeurs estimées par les splines pour chaque horizon standard 
  #             3 colonnes obligatoires: Soil.ID, Upper.Boundary, Lower.Boundary
  #             n colonnes : une colonne par propriété de sol 
  
  for (i in 1:length(champ.variables)){
    pred.i <- pred.spl.multiP2(data, std.depth, champ.id, champ.prof.sup, champ.prof.inf, champ.variables[i])
    pred.i <- pred.i[[1]]
    id.tot <- data.matrix(pred.i[,c(champ.id, champ.prof.sup, champ.prof.inf)],rownames.force=FALSE)
    id.tot <- as.matrix(apply(id.tot, 1, conc.champ))
    pred.i$id_tot <- id.tot 
    
    if (i == 1){
      res <-  pred.i
    } else {
      res <- merge(res, pred.i, by = "id_tot",all.x = F, all.y = F, suffixes = c("",".y"))
      res <- res[, -c(ncol(res)-1, ncol(res)-2, ncol(res)-3)]
    }
  }
  res <- res[,-1]
  return(res)
}

prep.hz <- function (table, champ.id, champ.prof.sup, champ.prof.inf, champ.variable, champ.n.horizon, champ.till.depth, champ.occup, liste.occup, default.till.width = 30, xd = 0.000001){
  # DESCRIPTION : Préparation des horizons pour spline : ajout de deux horizons d'épaisseur xd aux extrémités des horizons labourés
  # Cette fonction doit être appliquée à un seul profil
  # - labouré
  # - avec au moins deux horizons
  # - avec des valeurs de la propriété d'intérêt dans l'horizon labouré
  
  # On considère 5 organisations possibles  des horizons : voir fichier pdf associé pour la description des cas.
  #
  # N.B.: définition de la profondeur de labour :
  # - les horizons sont décrits (avec un horizon labouré) : on prend la profondeur max indiquée
  # - les horizons ne sont pas décrits et l'occupation du sol implique un horizon labourée : valeur par défaut indiquéee par l'argument "default.till.depth"
  # - les horizons ne sont pas décrit et l'occupation du sol est inconnue : on ne fait rien
  #  
  ## Inputs :
  # - table : soil horizon values (vector, mandatory)
  # - upper.boundary : upper limits of the described soil horizons (vector, mandatory)
  # - lower.boundary : upper limits of the described soil horizons (vector, mandatory)
  # - till.depth : profondeur de labour
  # - default.till.depth : profondeur de labour par défaut. Valeur par défaut : 30 cm.
  # - xd : épaisseur de l'horizon à ajouter aux extrémitées des horizons labourés
  
  ## Outputs : 
  # - data frame contenant la description des horizon
  
  table.r <- matrix(data = NA, nrow=nrow(table)*3,ncol=ncol(table))  #matrix a partir du data.frame originale
  table.r <- as.data.frame(table.r) #on retourne au data.frame
  colnames(table.r) <- colnames(table)
  
  # Définition de la profondeur inférieure du labour 
  if (!is.na(unique(table[,champ.till.depth]))){
    xlim <- unique(table[,champ.till.depth])
  } else {
    if (champ.occup %in% liste.occup){
 # ??? FIXME does not match the parameter list    xlim <- default.till.depth    
      xlim <- default.till.width   

    }
  }
  
  #Initialisation constantes et variables####################################################
  xdi  <- xd     #horizon infinitesimal 
  Pi   <- which( colnames(table) == champ.n.horizon)        #index propieté dans le tableau
  xres <- xlim     #profondeur residuale initiale = profondeur limite
  xsum <- 0        #sum in profondeur calcul moyenne entre 0 - xlim
  Yavg <- 0        #moyenne entre 0 - xlim et dans le boucle
  
  xi_r   <- 0      #memoire de valeurs si l'horizon est au milieu de la limite
  xii_r  <- 0
  chk_topsoil <- 0
  Yi_r   <- 0
  chk_prof  <- 1
  
  i <- 1
  n <- 0      #conteur ligne in ecriture
  hi <- 0
  chk_ir <- 0 #check horizon a cheval limite critique
  
  x_res <- FALSE  #boolean presence horizon residuel (apres xlim...)
  
  if ((all(table[,champ.occup] %in% liste.occup))) # on vérifie si le profil est labouré
    
  { 
    #-in case of labour---------------------------------------------------------
    
    #-------Start horizons LOOP-------------------------------------------------------
    for (i in 1:nrow(table))   # variable index i
    {
      #input dataset -> entre le single horizon ---------------------------------
      
      iDp  <- table[i,champ.id]
      iDh  <- table[i,champ.n.horizon]
      xi   <- table[i,champ.prof.sup]
      xii  <- table[i,champ.prof.inf]
      Occ  <- table[i,champ.occup]
      ptill <- table[i,champ.till.depth]
      
      Yi   <- table[i,champ.variable]
      
      if (xi < xlim)                   # 1) si le profondeur superieur de la couche est < xlim (e.g. 30cm)
      {
        iDh_1 <- iDh                     # le numero du premier horizon va en memoire
        
        if (xii > xlim)                  # 2) si le profondeur inferieure de la couche est > xlim (e.g. 30cm) -> on memorise la couche decoupé
        {                                   # 
          
          iDp_r  <- table[i,champ.id]  #id_profil
          iDh_r  <- table[i,champ.n.horizon]   #id_horizon
          chk_ir <- 1        #check decoupage -> (boolean) on memorise qu'on a decoup? l'horizon
          xi_r   <- xlim      #limite
          xii_r  <- xii  	  #fin de l'horizon
          Yi_r   <- Yi		    #proprieté principale 
          
          xii <- xlim         #decaler sur la limite -> le nouveau limite inferieure est pos? sur xlim
          x_res <- TRUE
        }
        #calcul de les moyennes
        Yavg <- (Yavg+(Yi*(xii-xi)/xlim))     #ponderation de l'horizon i
        xsum <- (xii-xi)+xsum	              #somme des horizons utilis?s en profondeur
        xres <- xres-(xii-xi)                #residual des des horizonse entre 0 - xlim
      }
      
      if ((xii >= xlim) || (x_res == TRUE))                     # 3) on a pass? l'horizon limite
      {
        if (xsum > 0)                                             #avec des valeurs en memoire
          #retablire la valeur moyenne s'il y a des horizons vides
        {
          Yavg=Yavg+(Yavg*(xlim/xsum))*xres/xlim
          
          # 4) ecrire les deux horizons limite 
          
          #Calcul la moyenne 0 - xlim on ecrit l'horizon 0 et l'horizon xlim
          #---------------------------------------------------------
          #ecriture horizon zero et xlimite (2 lignes)
          
          iDp <- table[i,champ.id]
          iDh  <- table[i,champ.n.horizon]
          xi   <- table[i,champ.prof.sup]
          xii  <- table[i,champ.prof.inf]
          Occ  <- table[i,champ.occup]
          ptill  <- table[i,champ.till.depth]
          
          hi <- 1
          n <- n + 1             #conteur ligne in écriture
          table.r[n,champ.id] <- iDp             #id_profil
          table.r[n,champ.n.horizon] <- 1 #iDh_1           #iDp #+ "-1"   #id_horizon <----------------------
          table.r[n,champ.prof.sup] <- 0.0           #prof_sup
          table.r[n,champ.prof.inf] <- xdi           #prof_inf
          table.r[n,champ.variable] <- Yavg         #valeur proprieté horizon i
          table.r[n,champ.occup] <- Occ
          table.r[n,champ.till.depth] <- ptill
          
          hi <- 2
          n <- n + 1             #conteur ligne in écriture
          table.r[n,champ.id] <- iDp             #id_profil
          table.r[n,champ.n.horizon] <- 1 #iDh_1           #iDp #+ "-1"   #id_horizon <-----------------------
          table.r[n,champ.prof.sup] <- xdi           #prof_sup
          table.r[n,champ.prof.inf] <- xlim - xdi    #prof_inf
          table.r[n,champ.variable] <- Yavg        #valeur proprieté horizon i
          table.r[n,champ.occup] <- Occ
          table.r[n,champ.till.depth] <- ptill
          
          hi <- 3
          n <- n + 1             #conteur ligne in écriture
          table.r[n,champ.id] <- iDp             #id_profil
          table.r[n,champ.n.horizon] <- 1 #iDh_1           #iDp #+ "-1"   #id_horizon <------------------------
          table.r[n,champ.prof.sup] <- xlim - xdi    #prof_sup
          table.r[n,champ.prof.inf] <- xlim          #prof_inf
          table.r[n,champ.variable] <- Yavg         #valeur proprieté horizon i
          table.r[n,champ.occup] <- Occ
          table.r[n,champ.till.depth] <- ptill
          
          chk_topsoil <- 1       
          
          Yavg <- 0
          
          xsum <- 0
          xres <- xlim
          #---------------------------------------------------<<<<<<<<<
          # 5) Si on a memoire de l'horizon residuel on l'écrit
          #---------------------------------------------------<<<<<<<<<
          if (chk_ir== 1)
          {
            if (chk_topsoil == 1)
            {
              
              hi <- hi+1
              n <- n + 1             #conteur ligne in ?criture
              table.r[n,champ.id] <- iDp_r           #id_profil
              table.r[n,champ.n.horizon] <- iDh_r          #iDp #+ "-1"   #id_horizon
              
              table.r[n,champ.prof.sup] <- xi_r          #prof_sup
              table.r[n,champ.prof.inf] <- xii_r         #prof_inf
              table.r[n,champ.variable] <- Yi_r         #valeur propriet? horizon i
              table.r[n,champ.occup] <- Occ
              table.r[n,champ.till.depth] <- ptill
              
              #vider la memoire de l'horizon residuel
              iDp_r  <- 0
              iDp_r  <- 0
              xi_r   <- 0
              xii_r  <- 0
              Yi_r   <- 0
              
              chk_ir <- 0
              chk_topsoil <- 0
              x_res <- FALSE
              
            }
          }
        }
        else
        {
          # 6) Ecrire les horizons de profondit?
          #---------------------------------------------------------
          hi <- hi+1
          n <- n + 1             #conteur ligne in ecriture
          
          table.r[n,champ.id] <- iDp           #id_profil
          table.r[n,champ.n.horizon] <- iDh           #iDp #+ "-1"   #id_horizon
          table.r[n,champ.prof.sup] <- xi            #prof_sup
          table.r[n,champ.prof.inf] <- xii           #prof_inf
          table.r[n,champ.variable] <- Yi           #valeur propriet? horizon i 
          table.r[n,champ.occup] <- Occ
          table.r[n,champ.till.depth] <- ptill
          
          chk_ir <- 0
          chk_topsoil <- 0
          
        }
      }  #-------------------fin if ((xii >= xlim) || (x_res == TRUE)) 
      
    }  #---------------------------------------occup sol
  } else {#------------------------------------CICLE FOR
    n <- n +  1
    table.r[n,] <- table[i,]
  }
  table.r <- table.r[is.finite(table.r[,champ.id]),]
  return(table.r)
}


# VALIDATION FUNCTION
# TODOS : add a test and an example for the multiple stratum case
#######################################################################
#' Préparation des horizons pour spline : ajout de deux horizons d'épaisseur xd aux
#' extrémités des horizons labourés
#' Cette fonction doit être appliquée à plusieurs profils
#' @param table : soil horizon values (vector, mandatory)
#' @param till.depth : profondeur de labour
#' @param default.till.width : profondeur de labour par défaut. Valeur par défaut : 30 cm.
#' @param xd : épaisseur de l'horizon à ajouter aux extrémitées des horizons labourés
#' @param champ.id  TODO
#' @param champ.prof.sup TODO
#' @param champ.prof.inf TODO
#' @param champ.variable TODO
#' @param champ.n.horizon TODO
#' @param champ.till.depth TODO
#' @param champ.occup TODO
#' @param liste.occup TODO
#' @return  data frame contenant la description des horizon
#' @keywords profile processing
prep.hz.multi <- function (table, champ.id, champ.prof.sup, champ.prof.inf, champ.variable, champ.n.horizon, champ.till.depth, champ.occup, liste.occup, default.till.width = 30, xd=0.000001){
  # DESCRIPTION : Préparation des horizons pour spline : ajout de deux horizons d'épaisseur xd aux extrémités des horizons labourés
  # Cette fonction doit être appliquée à plusieurs profils
  ## Inputs :
  # - table : soil horizon values (vector, mandatory)
  # - upper.boundary : upper limits of the described soil horizons (vector, mandatory)
  # - lower.boundary : upper limits of the described soil horizons (vector, mandatory)
  # - till.depth : profondeur de labour
  # - default.till.depth : profondeur de labour par défaut. Valeur par défaut : 30 cm.
  # - xd : épaisseur de l'horizon à ajouter aux extrémitées des horizons labourés
  ## Outputs : 
  # - data frame contenant la description des horizon
  
  liste <- unique(table[,champ.id])
  # Initialisation de la barre de progression  
  pb <- txtProgressBar(min = 0, max =length(liste) , style = 3)
  
  for (i in 1: length(liste)){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    table.i <- table[table[,champ.id] == liste[i],]
    res.i <- prep.hz(table.i, champ.id, champ.prof.sup, champ.prof.inf, champ.variable, champ.n.horizon,
                     champ.till.depth, champ.occup, liste.occup, default.till.width = 30, xd)
    
    if (i == 1){
      res <- res.i
    } else {
      res <- rbind(res, res.i)
    }                  
  }
  close(pb)
  return(res)
}

max.min.var <- function(variable) {
  ## DESCRIPTION : function to determine the max and min of the analyte for plotting
  ## INPUTS
  # - variable : vecteur contenant la valeur de la propriété de sol d'intérêt
  ## OUTPUTS : 
  # - retval : liste de deux éléments (deux vecteurs de longueur 1). Le premier élément contient la valeur maximale, le second la valeur minimale.
  
  plotmax<- 1.1*max(variable)    # max conc of analyte
  #plotspline<-1.1*max(tspline_a)   # max conc of spline
  #varspline<- c(plotmax, plotspline)  # vector of plotmax and plotspline
  #xmax<-max(varspline)         # max of x axis
  
  plotmin<- 0.4*min(variable)    # min conc of analyte
  #splinemin<-0.4*min(tspline_a)   # min conc of spline
  #vectmin<- c(plotmin, splinemin)  # vector of plotmax and plotspline
  #xmin<-min(vectmin)         # mix of x axis
  retval <- list(plotmax, plotmin)
  return(retval)
}

plot_raw <- function(prof.sup, prof.inf, variable, var.name){
  # DESCRIPTION : fonction permettant de représenter graphiquement les valeurs observées sur un profil
  # INPUTS
  # - prof.sup : vecteur contenant les limites supérieures des horizons
  # - prof.inf : vecteur contenant les limites inférieures des horizons  
  # - variable : vecteur contenant la valeur d'intérêt (valeur brute)
  # - var.name : nom de la variable d'intérêt (légende de l'axe des x sur le graphique en sortie)
  # OUTPUTS
  # - graphique
  
  t <- t(as.matrix(prof.sup)) # top of soil layer
  b <- t(as.matrix(prof.inf)) # bottom of soil layer
  c <- t(as.matrix(variable)) # conc of analyte
  #First -X axis -analyte concentration
  poly1 <- c(c[1]-c[1], c[1], c[1], c[1]-c[1]) #layer1
  poly2 <- c(c[2]-c[2], c[2], c[2], c[2]-c[2]) #layer2
  poly3 <- c(c[3]-c[3], c[3], c[3], c[3]-c[3]) #layer3
  poly4 <- c(c[4]-c[4], c[4], c[4], c[4]-c[4]) #layer4
  poly5 <- c(c[5]-c[5], c[5], c[5], c[5]-c[5]) #layer5
  poly6 <- c(c[6]-c[6], c[6], c[6], c[6]-c[6]) #layer6
  poly7 <- c(c[7]-c[7], c[7], c[7], c[7]-c[7]) #layer7
  poly8 <- c(c[8]-c[8], c[8], c[8], c[8]-c[8]) #layer8
  # Second-Y axis - Soil depth
  poly1y <- c(t[1], t[1], b[1], b[1]) #layer1
  poly2y <- c(t[2], t[2], b[2], b[2]) #layer2
  poly3y <- c(t[3], t[3], b[3], b[3]) #layer3
  poly4y <- c(t[4], t[4], b[4], b[4]) #layer4
  poly5y <- c(t[5], t[5], b[5], b[5]) #layer5
  poly6y <- c(t[6], t[6], b[6], b[6]) #layer6
  poly7y <- c(t[7], t[7], b[7], b[7]) #layer7
  poly8y <- c(t[8], t[8], b[8], b[8]) #layer8
  
  mm <- max.min.var(variable)
  raw.plot<-plot(poly1, poly1y,type="n",ylim=c(max(prof.inf),0),xlim=c(mm[[2]],mm[[1]]),ylab="Depth (cm)", xlab= var.name, lty=2, lwd=1, xaxs="i", col="black", font.lab=2,cex.lab=1.5)
  
  #Polygon
  polygon (poly1,poly1y, lty=1, lwd=3, border="black")
  polygon (poly2,poly2y, lty=1, lwd=3, border="black")
  polygon (poly3,poly3y, lty=1, lwd=3, border="black")
  polygon (poly4,poly4y, lty=1, lwd=3, border="black")
  polygon (poly5,poly5y, lty=1, lwd=3, border="black")
  polygon (poly6,poly6y, lty=1, lwd=3, border="black")
  polygon (poly7,poly7y, lty=1, lwd=3, border="black")
  polygon (poly8,poly8y, lty=1, lwd=3, border="black")
  
  return (raw.plot)
}

plot_raw.sp <- function(prof.sup, prof.inf, variable, spl.pc, var.name){
  # DESCRIPTION : fonction permettant de représenter graphiquement les valeurs observées sur un profil
  # INPUTS
  # - prof.sup : vecteur contenant les limites supérieures des horizons
  # - prof.inf : vecteur contenant les limites inférieures des horizons  
  # - variable : vecteur contenant la valeur d'intérêt (valeur brute)
  # - spl.pc : vecteur de longueur k+1 contenant les valeurs prédites par les splines de 0 à k cm
  #            spl.pc est obtenu par la fonction "predSpline", avec par exemple l'argument zs = seq(0,k,1)
  # - var.name : nom de la variable d'intérêt (légende de l'axe des x sur le graphique en sortie)
  # OUTPUTS
  # - graphique
  
  pspl <- spl.pc
  
  t <- t(as.matrix(prof.sup)) # top of soil layer
  b <- t(as.matrix(prof.inf)) # bottom of soil layer
  c <- t(as.matrix(variable)) # conc of analyte
  #First -X axis -analyte concentration
  poly1 <- c(c[1]-c[1], c[1], c[1], c[1]-c[1]) #layer1
  poly2 <- c(c[2]-c[2], c[2], c[2], c[2]-c[2]) #layer2
  poly3 <- c(c[3]-c[3], c[3], c[3], c[3]-c[3]) #layer3
  poly4 <- c(c[4]-c[4], c[4], c[4], c[4]-c[4]) #layer4
  poly5 <- c(c[5]-c[5], c[5], c[5], c[5]-c[5]) #layer5
  poly6 <- c(c[6]-c[6], c[6], c[6], c[6]-c[6]) #layer6
  poly7 <- c(c[7]-c[7], c[7], c[7], c[7]-c[7]) #layer7
  poly8 <- c(c[8]-c[8], c[8], c[8], c[8]-c[8]) #layer8
  # Second-Y axis - Soil depth
  poly1y <- c(t[1], t[1], b[1], b[1]) #layer1
  poly2y <- c(t[2], t[2], b[2], b[2]) #layer2
  poly3y <- c(t[3], t[3], b[3], b[3]) #layer3
  poly4y <- c(t[4], t[4], b[4], b[4]) #layer4
  poly5y <- c(t[5], t[5], b[5], b[5]) #layer5
  poly6y <- c(t[6], t[6], b[6], b[6]) #layer6
  poly7y <- c(t[7], t[7], b[7], b[7]) #layer7
  poly8y <- c(t[8], t[8], b[8], b[8]) #layer8
  
  mm <- max.min.var(variable)
  
  raw.plot<-plot(poly1, poly1y,type="n",ylim=c(max(prof.inf),0),xlim=c(mm[[2]],mm[[1]]),ylab="Depth (cm)", xlab= var.name, lty=2, lwd=1, xaxs="i", col="black", font.lab=2,cex.lab=1.5)
  
  #Polygon
  polygon (poly1,poly1y, lty=1, lwd=3, border="black")
  polygon (poly2,poly2y, lty=1, lwd=3, border="black")
  polygon (poly3,poly3y, lty=1, lwd=3, border="black")
  polygon (poly4,poly4y, lty=1, lwd=3, border="black")
  polygon (poly5,poly5y, lty=1, lwd=3, border="black")
  polygon (poly6,poly6y, lty=1, lwd=3, border="black")
  polygon (poly7,poly7y, lty=1, lwd=3, border="black")
  polygon (poly8,poly8y, lty=1, lwd=3, border="black")
  yas<- as.matrix(t(c(seq(0,(length(pspl)-1),by=1))))
  lines(pspl, yas, lty=1, lwd=1,col="black")
  
  return (raw.plot)
}

plotAll <- function(id, prof.sup, prof.inf, variable, d, spl.pc, spl.h, var.name){
  # DESCRIPTION : fonction permettant de représenter graphiquement les valeurs observées sur un profil, plus les valeurs prédites
  #               par le spline  (par centimètres et pour des horizons standards)
  # INPUTS
  # - id : identifiant du profil étudié
  # - prof.sup : vecteur contenant les limites supérieures des horizons
  # - prof.inf : vecteur contenant les limites inférieures des horizons  
  # - variable : vecteur contenant la valeur d'intérêt (valeur brute)
  # - d : vecteurs de longeur n+1, décrivant les limites supérieures et inférieures des n horizons standards
  # - spl.pc : vecteur de longueur k+1 contenant les valeurs prédites par les splines de 0 à k cm
  #            spl.pc est obtenu par la fonction "predSpline", avec par exemple l'argument zs = seq(0,k,1)
  # - spl.hz : vecteur de longeur n contenant les valeurs prédites par les splines pour chaque horizon standard
  #             spl.hz est obtenu par la fonction intgSpline
  # - var.name : nom de la variable d'intérêt (légende de l'axe des x sur le graphique en sortie)
  # OUTPUTS
  # - graphique
  
  d <- matrix(d,1)
  v_nyfit <- spl.h
  pspl <- spl.pc
  
  t <- t(as.matrix(prof.sup)) # top of soil layer
  b <- t(as.matrix(prof.inf)) # bottom of soil layer
  c <- t(as.matrix(variable)) # conc of analyte
  # First -X axis -analyte concentration
  poly1 <- c(c[1]-c[1], c[1], c[1], c[1]-c[1]) #layer1
  poly2 <- c(c[2]-c[2], c[2], c[2], c[2]-c[2]) #layer2
  poly3 <- c(c[3]-c[3], c[3], c[3], c[3]-c[3]) #layer3
  poly4 <- c(c[4]-c[4], c[4], c[4], c[4]-c[4]) #layer4
  poly5 <- c(c[5]-c[5], c[5], c[5], c[5]-c[5]) #layer5
  poly6 <- c(c[6]-c[6], c[6], c[6], c[6]-c[6]) #layer6
  poly7 <- c(c[7]-c[7], c[7], c[7], c[7]-c[7]) #layer7
  poly8 <- c(c[8]-c[8], c[8], c[8], c[8]-c[8]) #layer8
  # Second-Y axis - Soil depth
  poly1y <- c(t[1], t[1], b[1], b[1]) #layer1
  poly2y <- c(t[2], t[2], b[2], b[2]) #layer2
  poly3y <- c(t[3], t[3], b[3], b[3]) #layer3
  poly4y <- c(t[4], t[4], b[4], b[4]) #layer4
  poly5y <- c(t[5], t[5], b[5], b[5]) #layer5
  poly6y <- c(t[6], t[6], b[6], b[6]) #layer6
  poly7y <- c(t[7], t[7], b[7], b[7]) #layer7
  poly8y <- c(t[8], t[8], b[8], b[8]) #layer8
  
  mm <- max.min.var(variable)
  
  raw.plot <-plot(poly1, poly1y,type="n",ylim=c(max(d),0),xlim=c(mm[[2]],mm[[1]]),ylab="Depth (cm)", xlab= var.name,main= paste("Point",id, sep=" "), lty=2, lwd=1, xaxs="i", col="black", font.lab=2,cex.lab=1.5)
  
  p_depth <- max(d)
  mat_nyfit <- t(as.matrix(v_nyfit[1:(length(v_nyfit))]))
  m_nyfit <- t(mat_nyfit) # change number for each plot
  upper <- as.matrix(d[,1:length(m_nyfit)])
  lower <- as.matrix(d[,2:(length(m_nyfit)+1)])
  vari_t <- as.data.frame(cbind(upper,lower,m_nyfit))
  vari.mat <- matrix(NA,ncol=length(vari_t[1,]),nrow=length(vari_t[,1]))
  for (vt in 1:length(vari.mat[,1])) {
    vari.mat[vt,1]<- vari_t[vt,1]
    vari.mat[vt,3]<- vari_t[vt,3]
    if (p_depth > vari_t[vt,1] & p_depth < vari_t[vt,2])
    {vari.mat[vt,2]<- p_depth} else {vari.mat[vt,2]<- vari_t[vt,2]}}  
  vari_d <- as.data.frame(vari.mat)
  
  t1 <- t(vari_d[1]) # top of soil layer
  b1 <- t(vari_d[2]) # bottom of soil layer
  c1 <- t(vari_d[3]) # conc of analyte
  #First -X axis -analyte concentration
  poly11 <- c(c1[1]-c1[1], c1[1], c1[1], c1[1]-c1[1]) #layer1
  poly21 <- c(c1[2]-c1[2], c1[2], c1[2], c1[2]-c1[2]) #layer2
  poly31 <- c(c1[3]-c1[3], c1[3], c1[3], c1[3]-c1[3]) #layer3
  poly41 <- c(c1[4]-c1[4], c1[4], c1[4], c1[4]-c1[4]) #layer4
  poly51 <- c(c1[5]-c1[5], c1[5], c1[5], c1[5]-c1[5]) #layer5
  poly61 <- c(c1[6]-c1[6], c1[6], c1[6], c1[6]-c1[6]) #layer6
  poly71 <- c(c1[7]-c1[7], c1[7], c1[7], c1[7]-c1[7]) #layer7
  poly81 <- c(c1[8]-c1[8], c1[8], c1[8], c1[8]-c1[8]) #layer8
  # Second-Y axis - Soil depth
  poly1y1 <- c(t1[1], t1[1], b1[1], b1[1]) #layer1
  poly2y1 <- c(t1[2], t1[2], b1[2], b1[2]) #layer2
  poly3y1 <- c(t1[3], t1[3], b1[3], b1[3]) #layer3
  poly4y1 <- c(t1[4], t1[4], b1[4], b1[4]) #layer4
  poly5y1 <- c(t1[5], t1[5], b1[5], b1[5]) #layer5
  poly6y1 <- c(t1[6], t1[6], b1[6], b1[6]) #layer6
  poly7y1 <- c(t1[7], t1[7], b1[7], b1[7]) #layer7
  poly8y1 <- c(t1[8], t1[8], b1[8], b1[8]) #layer8
  
  polygon (poly11,poly1y1, lty=2, lwd=1, border="grey50", col="grey")
  polygon (poly21,poly2y1, lty=2, lwd=1, border="grey50", col="grey")
  polygon (poly31,poly3y1, lty=2, lwd=1, border="grey50", col="grey")
  polygon (poly41,poly4y1, lty=2, lwd=1, border="grey50", col="grey")
  polygon (poly51,poly5y1, lty=2, lwd=1, border="grey50", col="grey")
  polygon (poly61,poly6y1, lty=2, lwd=1, border="grey50", col="grey")
  polygon (poly71,poly7y1, lty=2, lwd=1, border="grey50", col="grey")
  polygon (poly81,poly8y1, lty=2, lwd=1, border="grey50", col="grey")
  
  #Polygon
  polygon (poly1,poly1y, lty=1, lwd=1.5, border="black")
  polygon (poly2,poly2y, lty=1, lwd=1.5, border="black")
  polygon (poly3,poly3y, lty=1, lwd=1.5, border="black")
  polygon (poly4,poly4y, lty=1, lwd=1.5, border="black")
  polygon (poly5,poly5y, lty=1, lwd=1.5, border="black")
  polygon (poly6,poly6y, lty=1, lwd=1.5, border="black")
  polygon (poly7,poly7y, lty=1, lwd=1.5, border="black")
  polygon (poly8,poly8y, lty=1, lwd=1.5, border="black")
  yas <- as.matrix(t(c(seq(0,(length(pspl)-1),by=1))))
  lines(pspl, yas, lty=1, lwd=2,col="red3")
  axis(2, at = NULL, labels = F, tick = T, col = "black")
  return (raw.plot)
}

plotSeq<- function (data, champ.id, champ.prof.sup, champ.prof.inf, champ.variable, d, pspl, strt ,fini, var.name, nom){
  # DESCRIPTION : Fonction permettant de représenter graphiquement les résultats des splines pour plusieurs profils 
  #               et d'exporter le résultat sous forme d'un pdf
  # INPUTS
  # - data : matrice contenant au moins 4 colonnes avec : les identifiants des profils, les limites supérieures des horizons,
  #         les limites inférieures des horizons, la valeur d'intérêt (valeur brute)
  # - nyfit : data.frame n ligne (pour n profils) et 9 colonnes contenant l'id du profil, le lambda utilisé pour le profil, les valeurs
  #           prédites par les splines pour chaque horizon standard et la profondeur maximale des horizons standards.
  #           nyfit est obtenu par la fonction "pred.spl.multiP1" (premier élement de la liste obtenu en sortie de "pred.spl.multiP1")
  # - spfit : data.frame de n colonnes pour n profils, et de x lignes contenant les valeurs prédites par les splines de 0 à x cm
  #          spfit est obtenu par la fonction pred.spl.multiP1" (second élement de la liste obtenu en sortie de "pred.spl.multiP1")
  # - strt : premier profil à représenter
  # - fini : dernier profil à représenter
  # - nom : nom du fichier pdf (ex. "splines.pdf")
  # OUTPUTS
  # - fichier pdf
  
  sp_dat <- split(as.data.frame(data), data[,1])
  
  pdf(nom)
  for (i in strt:fini) {           
    vari <- sp_dat[[i]]  
    spl.h <- pspl[[1]][i,]
    spl.h <- t(spl.h[-c(1,2,length(spl.h))])
    spl.pc <- pspl[[2]][,i]
    
    prof.sup <- vari[,champ.prof.sup]
    prof.inf <- vari[,champ.prof.inf]
    variable <- vari[,champ.variable]
    id <- unique(vari[,champ.id])
    
    plotAll(id, prof.sup, prof.inf, variable, d, spl.pc, spl.h, var.name)
    dev.next()
  }
  graphics.off()
}

conc.champ <- function(x){
  res <- paste(x,collapse=" ")
}


### non spline funciton but still usefull for integrating profiles 

### 
### orderedProfData : data.frame with at least a depth variable (real) and a property
### should be ordered as a function of hz number
### variable which name is passed with the propr argument
### the depth variable is indeed the with of the horizons. Their
### order indicate the depth
### prop : name of the property the function should return a value of
### z : depth of the prediction
predSoilProperty <- function(z, orderedProfData, prop = "clay"){
	h <- 1
	cumDepth <- 0
	while (z >= cumDepth){
		proValue <- orderedProfData[h, prop]
		cumDepth <- cumDepth + orderedProfData$width[h]
		if (h == dim(orderedProfData)[1]) return(proValue)
		h <- h + 1
		}
	return(proValue)
}
### testing
test <- data.frame(width = c(0, 12, 5, 4), clay = c(2, 30, 50, 20))
#expect_equal(predSoilProperty(200, test), 20)
#expect_equal(predSoilProperty(20, test), 20)
#expect_equal(predSoilProperty(17, test), 20)
#expect_equal(predSoilProperty(0, test), 30)

### Aplplique la fonction predSoilProperty sur un vecteur de profondeurs
### renvoie un vecteur de valeurs pour la propriété étudiée
### permet d'intégrer la proprété sur une profondeur donnée
predPropVector <- Vectorize(predSoilProperty, "z")
	#function(zs, profData, prop = "clay"){
	#sapply(zs, predSoilProperty, profData, prop)
	#toto <- Vectorize(predSoilProperty, "z")
#}
predPropVectorOld <- #Vectorize(predSoilProperty, "z")
	function(zs, profData, prop = "clay"){
	sapply(zs, predSoilProperty, profData, prop)
	#toto <- Vectorize(predSoilProperty, "z")
}

depth <- 18


#' Integre une propriete sur un profil continu
#'
#' Pour un profil sans horizons manquants et pour une propriété donnée
#'
#' @param z profondeur d'intégration
#' @param orderedProfData dataframe contenant deux champs : le champ *width* des
#' épaisseurs des horizons successifs et un champ correspondant à prop.
#' Les lignes doivent être ordonnées de façon à correspondre aux horizons depuis la surface jusqu'au fond. 
#' @param prop nom de la propriété lower.boundary : upper limits of the described soil horizons (vector, mandatory)
#' @param widthName : variable donnant l'épaisseur des horizons
#' @return une valeur correspondant à l'intégration (la somme) de la propriété sur
#' l'épaisseur souhaitée.
#' @keywords  profile depth soil property integration
#' @export 
#' @examples
#' test <- data.frame(width = c(0, 12, 5, 4), clay = c(2, 30, 50, 20), carbone = c(0, 1.7, 0.5, 0.3))
#' depth <- 18
#' integrateContinuousLayers(depth, test, prop = "clay", widthName = "width")/depth
#' integrateContinuousLayers(depth, test, prop = "carbone", widthName = "width")/depth


integrateContinuousLayers <- function(z, orderedProfData, prop = "clay", 
					widthName = "depth"){
	h <- 1
	cumDepth <- 0#orderedProfData$width[h]
	proValue <- 0
	while ((z >= cumDepth + orderedProfData[h, widthName]) &&
		h <= (dim(orderedProfData)[1] - 1)){
		currDepth <- orderedProfData[h, widthName]
		cumDepth <- cumDepth + currDepth
		proValue <- proValue + orderedProfData[h, prop] * currDepth
#		if (h == dim(orderedProfData)[1]) return(proValue)

		h <- h + 1
		#browser()
	}
	remainingWitdh <- z - cumDepth
	proValue <- proValue + orderedProfData[h, prop] * remainingWitdh
	
	return(proValue)
}



