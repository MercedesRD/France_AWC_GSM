stats <- function(x){
  # Calculation of summary statistics (min, max, percentiles, median, mean, sd, and variance)
  stats <- matrix(NA,1,10); colnames(stats) <- c("Min","D1","Q1","Median","Mean","Q3","D9","Max","sd","var")
  stats[1,"Min"] <- min(x,na.rm=TRUE)
  stats[1,"Max"] <- max(x,na.rm=TRUE)
  stats[1,"D1"] <- quantile(x, seq(0,1,0.1),na.rm=TRUE)[2]  
  stats[1,"D9"] <- quantile(x, seq(0,1,0.1),na.rm=TRUE)[10]  
  stats[1,"Q1"] <- quantile(x, seq(0,1,0.25),na.rm=TRUE)[2]  
  stats[1,"Q3"] <- quantile(x, seq(0,1,0.25),na.rm=TRUE)[4]
  stats[1,"Median"] <- median(x,na.rm=TRUE)  
  stats[1,"Mean"] <- mean(x,na.rm=TRUE)
  stats[1,"sd"] <- sd(x,na.rm=TRUE)
  stats[1,"var"] <- var(x,na.rm=TRUE)
  stats <- round(stats,2)
  return(stats)
}

extract.data <- function(variable,login, password) {
  # variable : données à extraire, 5 arguments possibles : 
  # "analyses" --> données de la table data.resultat_analyse (toutes les analyses)
  # "da" --> données de la table data.resultat_densite_apparente
  # "eg.mes" --> données de la table data.resultat_eg (éléments grossiers mesurés)
  # "eg.surf" --> données de la table data.profils (éléments grossiers de surface; %)
  # "eg.est" --> données de la table data.horizon (éléments grossiers estimés)
  # "granulo" --> données de la table data.resultat_analyse_granulo
  # "prof.sol" --> données de la table data.profils (profondeur du sol en cm, profondeur des racines et raison d'arrêt)
  # "prof.sol.hz" --> données de la table data.profils (profondeur du sol en cm) + données de la table data.horizon (nom des horizons) 
  # "racines" --> données de la table data.horizon (abondance des racines)
  # "drainage" --> données de la table data.profils (drainage naturel)
  # "hydro" --> données de la table data.horizon (nature des taches et abondance)
  # "couleur" --> données de la table data.horizon (couleur des horizons)
  # "prelevements" --> informations sur les prélévements
  # "analyses" --> tous les résultats d'analyse chimique
  # Analyseschimiques :
  # "c_n", ph_eau", "ca_ech", "mg_ech", "k_ech", "cec","calc_tot", "fe_lib", "fe_tot", "n_tot", "p_tot", "al_tot",
  # "al_lib", "al_ech", "na_ech", "mn_ech", "h_ech", "fe_ech", "p_ass"
  
  require(RODBC)
  
  ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)
  
  if (variable == "da"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT resultat_densite_apparente.id_resultat, horizon.id_profil, horizon.no_horizon, horizon.prof_sup_min,
                             horizon.prof_sup_max, horizon.prof_sup_moy, horizon.prof_inf_min, horizon.prof_inf_max, 
                             horizon.prof_inf_moy, prelevement.id_prelevement, prelevement.no_prelevement, prelevement.prof_sommet, 
                             prelevement.prof_base, methode_analyse_physique.no_methode, methode_analyse_physique.nom_methode, resultat_densite_apparente.valeur,
                             horizon.cpcs_nom, horizon.rp_2008_nom, horizon.rp_95_nom
                             FROM data.resultat_densite_apparente
                             LEFT JOIN meta.methode_analyse_physique USING (id_methode)
                             JOIN data.analyse USING (id_analyse)
                             JOIN data.prelevement USING (id_prelevement)
                             JOIN data.horizon ON horizon.id_profil = prelevement.id_profil AND horizon.no_horizon = prelevement.no_horizon
                             WHERE data.analyse.id_prelevement IS NOT NULL
                             UNION 
                             SELECT resultat_densite_apparente.id_resultat, horizon.id_profil, horizon.no_horizon, horizon.prof_sup_min, horizon.prof_sup_max, 
                             horizon.prof_sup_moy, horizon.prof_inf_min, horizon.prof_inf_max, horizon.prof_inf_moy, NULL::bigint AS id_prelevemen,
                             NULL::text AS no_prelevement, NULL::integer AS prof_sommet, NULL::integer AS prof_base, methode_analyse_physique.no_methode,
                             methode_analyse_physique.nom_methode, resultat_densite_apparente.valeur,horizon.cpcs_nom, horizon.rp_2008_nom, horizon.rp_95_nom
                             FROM data.resultat_densite_apparente
                             LEFT JOIN meta.methode_analyse_physique USING (id_methode)
                             JOIN data.analyse USING (id_analyse)
                             JOIN data.horizon ON horizon.id_profil = data.analyse.id_profil AND data.horizon.no_horizon = data.analyse.no_horizon
                             WHERE data.analyse.id_prelevement IS NULL")    
  }
  
  if (variable == "eg.mes"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT prelevement.id_profil, prelevement.no_horizon, prelevement.id_prelevement, resultat_eg.prc_mass_tot, resultat_eg.masse, resultat_eg.volume, prof_sommet, prof_base
                     FROM data.resultat_eg
                     JOIN data.analyse USING (id_analyse)
                     JOIN data.prelevement ON prelevement.id_prelevement = data.analyse.id_prelevement
                     JOIN data.profil ON profil.id_profil = prelevement.id_profil
                     WHERE resultat_eg.prc_mass_tot IS NOT NULL;")
    
  }
  
  if (variable == "eg.est"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_profil, no_horizon, abondance_eg, abondance_eg_prin, abondance_eg_sec, nom_eg1_h, nom_eg2_h, taille_eg1_h, taille_eg2_h                     
                      FROM data.horizon
                      WHERE abondance_eg IS NOT NULL;")
  }
  
  if (variable == "eg.surf"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_profil, eg_surface
                     FROM data.profil
                     WHERE eg_surface IS NOT NULL;")
  }
  
  
  if (variable == "granulo"){
    sqlQuery(ds3,"SET client_encoding TO \'Latin1\'") 
    data <- sqlQuery(ds3,
                     "SELECT COALESCE(f3.id_prelevement, f2.id_prelevement, f1.id_prelevement) AS id_prelevement, f1.argile, f2.limon, f3.sable
                    FROM (SELECT r1.id_prelevement, avg(r1.valeur) AS sable
                     FROM (SELECT pr.id_prelevement, e.id_ens_resultat_analyse_granulo, 50 AS seuil_inferieur, 2000 AS seuil_superieur, 'µm' AS unite_seuil, sum(r.valeur) AS valeur
                     FROM data.ens_resultat_analyse_granulo e
                     JOIN data.analyse USING (id_analyse)
                     JOIN data.prelevement pr USING (id_prelevement)
                     JOIN data.profil p ON p.id_profil = pr.id_profil
                     JOIN data.resultat_analyse_granulo r USING (id_ens_resultat_analyse_granulo)
                     JOIN med.granulo_entre_990_et_1100 USING (id_ens_resultat_analyse_granulo)
                     LEFT JOIN data.preparation_granulo USING (id_ens_resultat_analyse_granulo)
                     WHERE (COALESCE(preparation_granulo.pretrait_d2, preparation_granulo.pretrait_d3) IS NULL OR (COALESCE(preparation_granulo.pretrait_d2, preparation_granulo.pretrait_d3) <> ALL (ARRAY['1'::text, '3'::text, '6'::text, '8'::text]))) AND r.seuil_inferieur >= 50::numeric AND r.seuil_inferieur < 2000::numeric AND r.seuil_superieur <= 2000::numeric AND r.seuil_superieur > 50::numeric AND r.valeur >= 0::numeric
                     GROUP BY pr.id_prelevement, e.id_ens_resultat_analyse_granulo) r1
                     GROUP BY r1.id_prelevement, r1.seuil_inferieur, r1.seuil_superieur, r1.unite_seuil::text) f3
                     FULL JOIN ( SELECT r1.id_prelevement, avg(r1.valeur) AS limon
                     FROM (SELECT pr.id_prelevement, e.id_ens_resultat_analyse_granulo, 2 AS seuil_inferieur, 50 AS seuil_superieur, 'µm' AS unite_seuil, sum(r.valeur) AS valeur
                     FROM data.ens_resultat_analyse_granulo e
                     JOIN data.analyse USING (id_analyse)
                     JOIN data.prelevement pr USING (id_prelevement)
                     JOIN data.profil p ON p.id_profil = pr.id_profil
                     JOIN data.resultat_analyse_granulo r USING (id_ens_resultat_analyse_granulo)
                     JOIN med.granulo_entre_990_et_1100 USING (id_ens_resultat_analyse_granulo)
                     LEFT JOIN data.preparation_granulo USING (id_ens_resultat_analyse_granulo)
                     WHERE (COALESCE(preparation_granulo.pretrait_d2, preparation_granulo.pretrait_d3) IS NULL OR (COALESCE(preparation_granulo.pretrait_d2, preparation_granulo.pretrait_d3) <> ALL (ARRAY['1'::text, '3'::text, '6'::text, '8'::text]))) AND r.seuil_inferieur >= 2::numeric AND r.seuil_inferieur < 50::numeric AND r.seuil_superieur <= 50::numeric AND r.seuil_superieur > 2::numeric AND r.valeur >= 0::numeric
                     GROUP BY pr.id_prelevement, e.id_ens_resultat_analyse_granulo) r1
                     GROUP BY r1.id_prelevement, r1.seuil_inferieur, r1.seuil_superieur, r1.unite_seuil::text) f2 ON f3.id_prelevement = f2.id_prelevement
                     FULL JOIN ( SELECT r1.id_prelevement, avg(r1.valeur) AS argile
                     FROM (SELECT pr.id_prelevement, e.id_ens_resultat_analyse_granulo, 0 AS seuil_inferieur, 2 AS seuil_superieur, 'µm' AS unite_seuil, sum(r.valeur) AS valeur
                     FROM data.ens_resultat_analyse_granulo e
                     JOIN data.analyse USING (id_analyse)
                     JOIN data.prelevement pr USING (id_prelevement)
                     JOIN data.profil p ON p.id_profil = pr.id_profil
                     JOIN data.resultat_analyse_granulo r USING (id_ens_resultat_analyse_granulo)
                     JOIN med.granulo_entre_990_et_1100 USING (id_ens_resultat_analyse_granulo)
                     LEFT JOIN data.preparation_granulo USING (id_ens_resultat_analyse_granulo)
                     WHERE (COALESCE(preparation_granulo.pretrait_d2, preparation_granulo.pretrait_d3) IS NULL OR (COALESCE(preparation_granulo.pretrait_d2, preparation_granulo.pretrait_d3) <> ALL (ARRAY['1'::text, '3'::text, '6'::text, '8'::text]))) AND r.seuil_inferieur >= 0::numeric AND r.seuil_inferieur < 2::numeric AND r.seuil_superieur <= 2::numeric AND r.seuil_superieur > 0::numeric AND r.valeur >= 0::numeric
                     GROUP BY pr.id_prelevement, e.id_ens_resultat_analyse_granulo) r1
                     GROUP BY r1.id_prelevement, r1.seuil_inferieur, r1.seuil_superieur, r1.unite_seuil::text) f1 ON f2.id_prelevement = f1.id_prelevement")
  
  # ajout de conditions sur les départements :  
  #  "(p.no_dep = ANY (ARRAY['45'::text, '18'::text, '36'::text, '37'::text, '41'::text, '28'::text])) "
  #  ajout de conditions sur les types de profils : 
  #  "(p.type_prof = ANY (ARRAY['2'::text, '3'::text, '5'::text]))"
    
  }
  
  if (variable == "carbone" | variable =="mat_org"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,paste("
         SELECT resultat_analyse.id_resultat, horizon.id_profil, horizon.no_horizon, horizon.prof_sup_min,
                horizon.prof_sup_max, horizon.prof_sup_moy, horizon.prof_inf_min, horizon.prof_inf_max, 
                horizon.prof_inf_moy, prelevement.id_prelevement, prelevement.no_prelevement, prelevement.prof_sommet, 
                prelevement.prof_base, methode_analyse.no_methode, methode_analyse.nom_methode, resultat_analyse.valeur,
                horizon.cpcs_nom, horizon.rp_2008_nom, horizon.rp_95_nom
         FROM data.resultat_analyse
         LEFT JOIN meta.methode_analyse USING (id_methode)
         JOIN data.analyse USING (id_analyse)
         JOIN data.prelevement USING (id_prelevement)
         JOIN data.horizon ON horizon.id_profil = prelevement.id_profil AND horizon.no_horizon = prelevement.no_horizon
         WHERE resultat_analyse.determination = '",variable,"'::text AND data.analyse.id_prelevement IS NOT NULL
  UNION 
        SELECT resultat_analyse.id_resultat, horizon.id_profil, horizon.no_horizon, horizon.prof_sup_min, horizon.prof_sup_max, 
                horizon.prof_sup_moy, horizon.prof_inf_min, horizon.prof_inf_max, horizon.prof_inf_moy, NULL::bigint AS id_prelevemen,
                NULL::text AS no_prelevement, NULL::integer AS prof_sommet, NULL::integer AS prof_base, methode_analyse.no_methode,
                methode_analyse.nom_methode, resultat_analyse.valeur,horizon.cpcs_nom, horizon.rp_2008_nom, horizon.rp_95_nom
        FROM data.resultat_analyse
        LEFT JOIN meta.methode_analyse USING (id_methode)
        JOIN data.analyse USING (id_analyse)
        JOIN data.horizon ON horizon.id_profil = data.analyse.id_profil AND data.horizon.no_horizon = data.analyse.no_horizon
        WHERE data.resultat_analyse.determination = '",variable,"'::text AND data.analyse.id_prelevement IS NULL;",sep=""))
    
    # Identification des données igcs (type de profil)
    # 1 : Profil fictif : profil théorique non réel décrit par les modales de U.T.S. ; sert à caractériser de façon théorique l’U.T.S.
    # 2 : Profil vrai : profil de sol réel (fosse, coupe, ...)
    # 3 : Sondage : sondage effectué à la tarière
    # 4 : Analyse agronomique : analyse des horizons de surface (sans véritable description du sol en profondeur)
    # 5 : Profil composite
    
    prof.igcs <- sqlQuery(ds3,"
         SELECT id_profil, type_prof
         FROM data.profil") 
        
    # Identification des données rmqs (type de profil)
    # type de profil :
    # B : composite Biosoil
    # C : Composite
    # D : prélèvement volumétrique
    # F : Fosse
    # M : 2ème composite
    # N : 3ème composite
    # O : 4ème composite
    # S : sondage
    # T : théorique
    # X : 2ème sondage
    # Y : 2ème fosse
    
    prof.rmqs <- sqlQuery(ds3,"
         SELECT id_profil, type_profil_rmqs
         FROM data_rmqs.l_profil_intervention") 

    programme <- merge(prof.igcs, prof.rmqs, by = "id_profil", all.x = T, all.y = T)
    programme$type_prof <- as.character(programme$type_prof)
    programme$type_profil_rmqs <- as.character(programme$type_profil_rmqs)
    programme$id_profil <- as.character(programme$id_profil)
    programme$type <- programme$type_prof 
    programme$type[is.na(programme$type)] <- programme$type_profil_rmqs[is.na(programme$type)]
    programme$type_prof <- as.factor(programme$type_prof)
    programme$type_profil_rmqs <- as.factor(programme$type_profil_rmqs)
    programme$type <- as.factor(programme$type)
    data <- merge(data, programme, by = "id_profil", all.x = T, all.y = F)
  }
  
  if (  variable == "al_ech" | variable == "al_lib" | variable == "al_tot" | variable == "c_n" | variable == "ca_ech" | variable == "ca_tot" |variable == "calc_tot" |
        variable == "cd_ext" | variable == "cd_tot" | variable == "cec" | variable == "co_ext" | variable == "co_tot" | variable == "cr_ext" | variable == "cr_tot" | 
        variable == "cu_ext" | variable == "cu_tot" | variable == "h_ech" | variable == "fe_ech" | variable == "fe_ext" |variable == "fe_lib" | variable == "fe_tot" | 
        variable == "k_ech" | variable == "k_tot" | variable == "mg_ech" | variable == "mg_tot" | variable == "mn_ech" |  variable == "mn_tot" | variable == "n_tot" |
         variable == "ni_ext" |variable == "ni_tot" |variable == "na_ech" | variable == "na_tot" |variable == "p_tot"|variable == "p_ass"| 
        variable == "pb_ext" |variable == "pb_tot" |  variable == "ph_eau" |variable == "teneur_eau_res" | variable == "zn_ext" |variable == "zn_tot" |
          # analyses physiques
        variable == "cond_elec" ){ 
  sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
  data <- sqlQuery(ds3,paste("
        SELECT data.analyse.id_profil, data.analyse.no_horizon, id_analyse, id_prelevement, prof_sommet, prof_base, id_resultat, id_methode, no_methode, nom_methode,  unite, determination, valeur
        FROM data.resultat_analyse
        INNER JOIN data.analyse USING(id_analyse)
        INNER JOIN meta.methode_analyse USING(id_methode)
        INNER JOIN data.prelevement USING(id_prelevement)
        Where determination = '",variable,"'",sep=""))
  }
  
  if (variable == "analyses"){ 
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
        "SELECT data.analyse.id_profil, data.analyse.no_horizon, id_analyse, id_prelevement, prof_sommet, prof_base, id_resultat, id_methode, no_methode, nom_methode,  unite, determination, valeur
        FROM data.resultat_analyse
        INNER JOIN data.analyse USING(id_analyse)
        INNER JOIN meta.methode_analyse USING(id_methode)
        INNER JOIN data.prelevement USING(id_prelevement)
        Where determination IN ('al_ech', 'al_lib', 'al_tot','c_n','ca_ech','ca_tot','calc_act','calc_tot','cd_ext','cd_tot','cec','co_ext','co_tot','cr_ext','cr_tot',
                                'cu_ext','cu_tot','h_ech','fe_ech','fe_ext','fe_lib','fe_tot',
                                'k_ech','k_tot','mg_ech','mg_tot','mn_ech','mn_tot','n_tot',
                                'ni_ext','ni_tot','na_ech','na_tot','p_tot','p_ass',
                                'pb_ext','pb_tot','ph_eau','teneur_eau_res','zn_ext','zn_tot',
                                'cond_elec' )")
  }
  
  
  if (variable == "prof.sol"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_profil, prof_sol_p, abond_rac_p, distrib_rac_p, arret
                     FROM data.profil
                     WHERE prof_sol_p >= 0 AND prof_sol_p IS NOT null;")
  }
    if (variable == "prof.sol.hz"){
      sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
      data <- sqlQuery(ds3,
                       "SELECT id_profil, no_horizon, data.horizon.cpcs_nom, data.horizon.rp_95_nom, data.horizon.rp_2008_nom, 
                        prof_sup_moy, prof_inf_moy, prof_sol_p
                        FROM data.horizon
                        INNER JOIN data.profil USING(id_profil);")
      
  }
  if (variable == "racines"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_profil, no_horizon, prof_sup_moy, prof_inf_moy, abond_rac_h
                     FROM data.horizon;")
  }
  if (variable == "drainage"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_profil, drai_nat_p
                     FROM data.profil
                     WHERE drai_nat_p IS NOT null;")
  }
  if (variable == "hydro"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_profil, no_horizon, prof_sup_moy, prof_inf_moy, nat_tach_1, nat_tach_2, nat_tach_3, abond_tach_1,abond_tach_2,abond_tach_3
                     FROM data.horizon;")
  }
  if (variable == "couleur"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_profil, no_horizon, prof_sup_moy, prof_inf_moy, coul_h
                     FROM data.horizon
                     WHERE coul_h IS NOT null;")
  }
  if (variable == "prelevements"){
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    data <- sqlQuery(ds3,
                     "SELECT id_prelevement, no_prelevement, prof_sommet, prof_base, id_profil, no_horizon 
                     FROM data.prelevement;")
  }
  
  
  if (variable == "horizons"){
      sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
      data <- sqlQuery(ds3,
                       "SELECT id_profil, no_horizon,cpcs_nom, rp_95_nom, rp_2008_nom, prof_sup_moy, prof_inf_moy
                     FROM data.horizon;")
                       
                       # Identification des données igcs (type de profil)
                       # 1 : Profil fictif : profil théorique non réel décrit par les modales de U.T.S. ; sert à caractériser de façon théorique l’U.T.S.
                       # 2 : Profil vrai : profil de sol réel (fosse, coupe, ...)
                       # 3 : Sondage : sondage effectué à la tarière
                       # 4 : Analyse agronomique : analyse des horizons de surface (sans véritable description du sol en profondeur)
                       # 5 : Profil composite
                       
     prof.igcs <- sqlQuery(ds3,"SELECT id_profil, type_prof
                           FROM data.profil") 
                       
                       # Identification des données rmqs (type de profil)
                       # type de profil :
                       # B : composite Biosoil
                       # C : Composite
                       # D : prélèvement volumétrique
                       # F : Fosse
                       # M : 2ème composite
                       # N : 3ème composite
                       # O : 4ème composite
                       # S : sondage
                       # T : théorique
                       # X : 2ème sondage
                       # Y : 2ème fosse
                       
    prof.rmqs <- sqlQuery(ds3,"SELECT id_profil, type_profil_rmqs
                          FROM data_rmqs.l_profil_intervention") 
                       
                       programme <- merge(prof.igcs, prof.rmqs, by = "id_profil", all.x = T, all.y = T)
                       programme$type_prof <- as.character(programme$type_prof)
                       programme$type_profil_rmqs <- as.character(programme$type_profil_rmqs)
                       programme$id_profil <- as.character(programme$id_profil)
                       programme$type <- programme$type_prof 
                       programme$type[is.na(programme$type)] <- programme$type_profil_rmqs[is.na(programme$type)]
                       programme$type_prof <- as.factor(programme$type_prof)
                       programme$type_profil_rmqs <- as.factor(programme$type_profil_rmqs)
                       programme$type <- as.factor(programme$type)
                       data <- merge(data, programme, by = "id_profil", all.x = T, all.y = F)
                         
                       
  }
  
  
  
  
  odbcCloseAll()
  return(data)   
}

lim.hz.neg <- function(table,nom){
  # N.B. : function pour supprimer les horizons avec des limites négatives pour prof_sommet et prof_base
  table.prof <- table[,c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")]
  # supression des horizons avec limite négative
  liste.del <- ((!is.na(table.prof[,"prof_sommet"]) & table.prof[,"prof_sommet"] < 0) | 
               (!is.na(table.prof[,"prof_base"]) & table.prof[,"prof_base"] < 0) | 
               (!is.na(table.prof[,"prof_sup_moy"]) & table.prof[,"prof_sup_moy"] < 0) |
               (!is.na(table.prof[,"prof_inf_moy"]) & table.prof[,"prof_inf_moy"] < 0) |
               (is.na(table.prof[,"prof_sommet"]) & !is.na(table.prof[,"prof_sup_moy"])  & table.prof[,"prof_sup_moy"] < 0) | 
               (is.na(table.prof[,"prof_base"]) & !is.na(table.prof[,"prof_inf_moy"]) & table.prof[,"prof_inf_moy"] < 0))
  prof.del <- subset(table, liste.del)
  write.table(prof.del,paste("lim_hz_neg_",nom,".csv",sep=""),sep=";")
  table <- subset(table, !(liste.del))
  return(table)
}

del.hz.hum <- function(table,nom.hz.hum,nom.table){
  # trois champs possibles pour les noms d'horizons : cpcs_nom, rp_2008_nom, rp_95_nom
  # horizons humus existants : "A00","A0F","A0H","OH","OF","OL"
  liste.hum <- vector()
  for (i in 1:length(nom.hz.hum))
  {
    hum <- grep(nom.hz.hum [i],table$cpcs_nom) 
    if (length(hum) != 0 ) liste.hum  <- c(liste.hum,hum)
    hum <- grep(nom.hz.hum [i],table$rp_2008_nom)
    if (length(hum) != 0 ) liste.hum  <- c(liste.hum,hum)
    hum <- grep(nom.hz.hum [i],table$rp_95_nom)
    if (length(hum) != 0 ) liste.hum  <- c(liste.hum,hum)
  }
  if (length(liste.hum) > 0){
  table.hum  <- table[liste.hum,]
  write.table(table.hum ,paste("del_hz_hum_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
  table.non.hum <- table[-liste.hum,]
  } else {
    print("pas d'horizon humiques")
    table.non.hum <- table
  }
  return(table.non.hum)
}

del.decal <- function(table, nom.table){
  table.del <- subset(table, prof_sup_moy > prof_sommet | prof_inf_moy < prof_base)
  table.cor <- subset(table, !(prof_sup_moy > prof_sommet | prof_inf_moy < prof_base))
  write.table(table.del ,paste("del_decal_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)  
  return(table.cor)
}

conc.champ <- function(x){
  res <- paste(x,collapse=" ")
}

del.duplicated1 <- function(table, champs, nom.table){ 
  # Permet de supprimer des lignes dont plusieurs champs sont dupliqués. Par exemple : mêmes limites de prelevement et même valeur mesurée
  champ.dupl <- data.matrix(table[,champs],rownames.force=FALSE)
  champ.dupl <- as.matrix(apply(champ.dupl, 1, conc.champ))
  liste.champ.dupl <- duplicated(champ.dupl)
  if (length(liste.champ.dupl[liste.champ.dupl == TRUE]) == 0) {
    print("pas de doublons")} else{
      prof.dupl <- table[liste.champ.dupl,]
      write.table(prof.dupl ,paste("del_duplicated1_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
      
      while(length(liste.champ.dupl[liste.champ.dupl == TRUE]) != 0){
        champ.dupl <- data.matrix(table[,champs],rownames.force=FALSE)
        champ.dupl <- as.matrix(apply(champ.dupl, 1, conc.champ))
        liste.champ.dupl <- duplicated(champ.dupl)
        table <- subset(table, !liste.champ.dupl)
      }
    }
  return(table)
}

del.duplicated2 <- function(table, nom.table, option, seuil = 1.5){ 
  # Permet de traiter des lignes donct les limites de prelevement et les méthodes d'analyses sont les mêmes, mais qui ont des valeurs mesurées différentes.
  # Ex : présence d'éléments particuliers dans une matrice (gloss, etc) : analyses faites sur la matrice et l'élément particulier.
  # Option :
  # del : supression de toutes les lignes dupliquées (pas moyen de savoir laquelle est la plus représentatiive de l'horizon)
  # mean : moyenne des valeurs
  # delif : supression des lignes dupliquées si la différence entre les mesures est > au seuil
  
  champs <- c("id_profil","prof_sommet","prof_base","no_methode")
  champ.dupl <- data.frame(table[,champs])
  champ.dupl <- as.matrix(apply(champ.dupl, 1, conc.champ))
  liste.champ.dupl <- duplicated(champ.dupl)
  if (length(liste.champ.dupl[liste.champ.dupl == TRUE]) == 0)
    {
    print("pas de doublons")} else {
      write.table(champ.dupl ,paste("del_duplicated2_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
      
      if (option == "del"){
        table <- subset(table, !liste.champ.dupl) 
      }
      
      if (option == "mean"){
        prof.dupl <- table[liste.champ.dupl,] 
        champs <- c("id_profil","no_horizon")
        id.tot <- data.matrix(prof.dupl[,champs],rownames.force=FALSE)
        id.tot <- as.matrix(apply(id.tot, 1, conc.champ))
        prof.dupl <- subset(prof.dupl, !duplicated(id.tot))
        
        for(i in 1:nrow(prof.dupl)) {
          hz.i <- prof.dupl[i,]
          liste.prof <- subset(table, id_profil == hz.i$id_profil &  prof_sommet == hz.i$prof_sommet & prof_base == hz.i$prof_base & no_methode == hz.i$no_methode)
          hz.i$valeur <- mean(liste.prof$valeur)
          table <- subset(table, !(id_profil == hz.i$id_profil & no_horizon == hz.i$no_horizon & no_methode == hz.i$no_methode))
          table <- rbind(table, hz.i)
        }
      }
      
      if (option == "delif"){
        prof.dupl <- table[liste.champ.dupl,]
        champs <- c("id_profil","no_horizon")
        id.tot <- data.matrix(prof.dupl[,champs],rownames.force=FALSE)
        id.tot <- as.matrix(apply(id.tot, 1, conc.champ))
        prof.dupl <- subset(prof.dupl, !duplicated(id.tot))
        
        for(i in 1:nrow(prof.dupl)){
          #print(i)
          hz.i <- prof.dupl[i,]
          liste.prof <- subset(table, id_profil == hz.i$id_profil &  prof_sommet == hz.i$prof_sommet & prof_base == hz.i$prof_base & no_methode == hz.i$no_methode)
          # liste.prof <- subset(table, id_profil == hz.i$id_profil & (no_horizon == hz.i$no_horizon | (prof_sommet == hz.i$prof_sommet & prof_base == hz.i$prof_base)) & no_methode == hz.i$no_methode)
          if(length(unique(liste.prof$id_prelevement)) != 1) {#print(i)
            min.i <- min(liste.prof$valeur)
            max.i <- max(liste.prof$valeur)
          
          if ((max.i - min.i) <= seuil){
            hz.i$valeur <- mean(liste.prof$valeur)
            table <- subset(table, !(id_profil == hz.i$id_profil & no_horizon == hz.i$no_horizon & no_methode == hz.i$no_methode))
            table <- rbind(table, hz.i) 
          } else {
            table <- subset(table, !(id_profil == hz.i$id_profil & no_horizon == hz.i$no_horizon & no_methode == hz.i$no_methode))
          }
        }
      }
      
    }
}
  return(table)
}

del.duplicated3 <- function(table, nom.table, seuil, champ.id){ 
  # permet de repérer des horizons avec les mêmes limites d'horizons, des valeurs différentes (par ex. carbone) et des méthodes d'analyses différentes
  # la méthode retenue est alors la plus couremment utilisée
  #champs <- c("id_profil","prof_sommet","prof_base")
  #champ.dupl <- data.matrix(table[,champs],rownames.force=FALSE)
  #champ.dupl <- as.matrix(apply(champ.dupl, 1, conc.champ))
  liste.champ.dupl <- duplicated(table[,champ.id])
  liste.champ.dupl <- table[liste.champ.dupl,champ.id]
  liste.champ.dupl <- unique(liste.champ.dupl)
  
  
  if (length(liste.champ.dupl) == 0){
    print("pas de doublons")} else {
      while (length(liste.champ.dupl) != 0){
        pb <- txtProgressBar(min = 0, max = length(liste.champ.dupl) , style = 3)
        
        prof.dupl <- table[liste.champ.dupl,]
        write.table(prof.dupl ,paste("del_duplicated3_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        
        methode <- orderBy(~x,aggregate(table$valeur, by = list(table$nom_methode), length)); colnames(methode) <- c("nom_methode","n")
        methode$n[methode$nom_methode == "Non connue"] <- 0
        for(i in 1:length(liste.champ.dupl))
        {
          Sys.sleep(0.1)
          setTxtProgressBar(pb, i)
          
          liste.prof <- table[table[,champ.id] == liste.champ.dupl[i],]
          #liste.prof <- subset(table, id_tot == liste.champ.dupl[i])
          methode.i <- unique(liste.prof$nom_methode)
          if (length(methode.i) > 1) {
            methode.ii <- subset(methode, nom_methode %in% methode.i)
            methode.ii <- subset(methode.ii, n == max(n))
            id <- subset(liste.prof, nom_methode != methode.ii$nom_methode)$id_resultat
            table <- subset(table, id_resultat %nin% id)
          } else {
            min.i <- min(liste.prof$valeur)
            max.i <- max(liste.prof$valeur)
            #table <- subset(table, id_tot != liste.champ.dupl[i])
            table  <- table[table[,champ.id] != liste.champ.dupl[i],]
            if (max.i - min.i <= seuil){
              liste.prof.i <- liste.prof[1,]
              liste.prof.i$valeur <- mean(liste.prof$valeur)
              table <- rbind(table, liste.prof.i)
              
            }
          }
        }
        close(pb)
        liste.champ.dupl <- duplicated(table[,champ.id])
        liste.champ.dupl <- table[liste.champ.dupl,champ.id]
        liste.champ.dupl <- unique(liste.champ.dupl)
       }         
      }
  table <- orderBy(~ id_profil + prof_sommet, table)
  return(table)
}
       
coherence.lim <- function(table,nom.table){
  # INPUT :

  table.prof <- table[,c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")]
  # supression des horizons ave limite négative
  liste.del <- (!is.na(table.prof[,"prof_sommet"]) & !is.na(table.prof[,"prof_base"]) & table.prof[,"prof_sommet"] >= table.prof[,"prof_base"]) | 
                ((is.na(table.prof[,"prof_base"]) | is.na(table.prof[,"prof_base"])) & (!is.na(table.prof[,"prof_sup_moy"]) & !is.na(table.prof[,"prof_inf_moy"])) 
                 & table.prof[,"prof_sup_moy"] >= table.prof[,"prof_inf_moy"]) 
  if (dim(subset(table, liste.del))[1] == 0){
    print("Pas de problèmes de cohérence de limite")} else {
      hz.probleme <-  subset(table, liste.del)
      write.table(hz.probleme ,paste("coherence_lim_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
      table <-  subset(table, !liste.del) 
    }
  return(table)
}

hz.superpose <- function(table, id, limup, limdown, champs.comp, id.comp){
  
  #permet de détecter des horizons qui se sperposent mais n'ont pas de limites communes
  liste <- unique(table[,id])
  pb <- txtProgressBar(min = 0, max = length(liste) , style = 3)
  sup.prof <- vector()
  sup.up <- vector()
  sup.dwn <- vector()
  for(i in 1:length(liste)){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    p.i <- table[table[,id]== liste[i],]
    p.i <- p.i[order(p.i[,limup],decreasing =F),]  
    if (nrow(p.i)>1){
      for (j in 2:((nrow(p.i)))){
        if (p.i[j,limup] < p.i[j-1,limdown]) {
          sup.prof <- c(sup.prof, p.i[j,id])
          sup.up <- c(sup.up , p.i[j,limup])
          sup.dwn <- c(sup.dwn, p.i[j,limdown])
        }
      }
    }
    }
  close(pb)
  if (length(sup.prof) > 0){
    sup <- cbind(sup.prof, sup.up, sup.dwn)
    colnames(sup) <- c(id, limup, limdown)
    sup <- as.data.frame(sup)
    sup <- add.id(sup, champs.comp, id.comp)
    table <- table[table[,id.comp] %nin% sup[,id.comp],]
    print("Présence d'horizons superposés")
    return(list(table,sup))
  } else {
    print("Pas d'horizons superposés")
    return(table)
  }
}

# functions lim.hz.NA.1, lim.hz.NA.2 and lim.hz.NA.3 : used to correct the fields of horizon limits
lim.hz.NA.1 <- function(table,nom.table){
  # N.B. : fonction adaptée pour repérer les erreurs pour les horizons prof_base, prof_sommet, prof_inf_moy et prof_sup_moy
  # Function used to detect the errors for the fields prof_base, prof_sommet, prof_inf_moy and prof_sup_moy
  if (dim(subset(table, is.na(prof_sommet) | is.na(prof_base)))[1]==0){ 
    print("Pas de NA dans les profondeurs de prélévement")} else {
      
    ### Assigns the limits of the horizons to the prelevements or viceversa, when the other one is missing
        
      # 0.1.  cas où prof_sommet= NA et prof_base = NA, avec prof_sup_moy != NA et prof_inf_moy != NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- is.na(prof[,"prof_sommet"]) & is.na(prof[,"prof_base"]) & !is.na(prof[,"prof_sup_moy"]) & !is.na(prof[,"prof_inf_moy"])
      table$prof_base[liste.na] <- table$prof_inf_moy[liste.na]
      table$prof_sommet[liste.na] <- table$prof_sup_moy[liste.na]
    
      # 0.2.  cas où prof_sommet != NA et prof_sup_moy = NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- !is.na(prof[,"prof_sommet"]) & is.na(prof[,"prof_sup_moy"])
      table$prof_sup_moy[liste.na] <- table$prof_sommet[liste.na]
      
      # 0.3.  cas où prof_base != NA et prof_inf_moy = NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- !is.na(prof[,"prof_base"]) & is.na(prof[,"prof_inf_moy"])
      table$prof_inf_moy[liste.na] <- table$prof_base[liste.na]
      
      # 0.4.  cas où prof_sommet != NA et prof_base != NA, avec prof_sup_moy = NA et prof_inf_moy = NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- !is.na(prof[,"prof_sommet"]) & !is.na(prof[,"prof_base"]) & is.na(prof[,"prof_sup_moy"]) & is.na(prof[,"prof_inf_moy"])
      table$prof_inf_moy[liste.na] <- table$prof_base[liste.na]
      table$prof_sup_moy[liste.na] <- table$prof_sommet[liste.na]
      ### when this happens and there are multiple prelevements by horizon, 
      ### we will have several mesurements assign to a same horizon id, but with different horizon limits
      
      # 1. cas où prof_sommet= NA, avec prof_sup_moy != NA 
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- is.na(prof[,"prof_sommet"]) & !is.na(prof[,"prof_sup_moy"])
      mat.na <- subset(table, liste.na)
      if (nrow(mat.na) !=0){ 
        write.table(mat.na,paste("lim_hz_NA_1a_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        table$prof_sommet[liste.na] <- table$prof_sup_moy[liste.na] ### we assign the horizon upper limit to the prelevement
      }
      
      # 2. cas où prof_sommet= NA et prof_sup_moy = NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- is.na(prof[,"prof_sommet"]) & is.na(prof[,"prof_sup_moy"])
      mat.na <- subset(table, liste.na)
      if (nrow(mat.na)  !=0)
      { 
        write.table(mat.na,paste("lim_hz_NA_1b_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        for(i in 1:nrow(mat.na)) {
          hz.i <- mat.na$no_horizon[i]
          if (hz.i != 1) { table <- subset(table, id_tot1 != mat.na$id_tot1[i])} else {
            table$prof_sommet[table$id_tot1 == mat.na$id_tot1[i]] <- 0   ### assigns 0 to the top, surface of the soil
            table$prof_sup_moy[table$id_tot1 == mat.na$id_tot1[i]] <- 0}
        }
      }
      
      # 3. cas où prof_base = NA et prof_inf_moy != NA et prof_sommet < prof_inf_moy
      ### When the lower limit of the prelevement is missing, and if the upper limit of the prelevement is shallower (upper)
      ### than the lower limit of the horizon, it assign the horizon limit to the prelevement lower limit
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- is.na(prof[,"prof_base"]) & !is.na(prof[,"prof_inf_moy"]) & (prof[,"prof_sommet"] < prof[,"prof_inf_moy"])
      mat.na <- subset(table, liste.na)
      if (nrow(mat.na) !=0){ 
        write.table(mat.na,paste("lim_hz_NA_1c_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        table$prof_base[liste.na] <- table$prof_inf_moy[liste.na] 
      }
      
      # 4. cas où prof_base = NA et prof_inf_moy != NA et prof_sommet >= prof_inf_moy
      ### When the lower limit of the prelevement is missing, and if the upper limit of the prelevement is equal or deeper
      ### than the lower limit of the horizon, it eliminates this horizon.
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- is.na(prof[,"prof_base"]) & !is.na(prof[,"prof_inf_moy"]) & (prof[,"prof_sommet"] >= prof[,"prof_inf_moy"])
      mat.na <- subset(table, liste.na)
      if (nrow(mat.na) !=0){ 
        write.table(mat.na,paste("lim_hz_NA_1d_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        table <- subset(table, !liste.na)  
      }
      
      # 0.1.  cas où prof_sommet= NA et prof_base = NA, avec prof_sup_moy != NA et prof_inf_moy != NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- is.na(prof[,"prof_sommet"]) & is.na(prof[,"prof_base"]) & !is.na(prof[,"prof_sup_moy"]) & !is.na(prof[,"prof_inf_moy"])
      table$prof_base[liste.na] <- table$prof_inf_moy[liste.na]
      table$prof_sommet[liste.na] <- table$prof_sup_moy[liste.na]
      
      # 0.2.  cas où prof_sommet != NA et prof_base != NA, avec prof_sup_moy = NA et prof_inf_moy = NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- !is.na(prof[,"prof_sommet"]) & !is.na(prof[,"prof_base"]) & is.na(prof[,"prof_sup_moy"]) & is.na(prof[,"prof_inf_moy"])
      table$prof_inf_moy[liste.na] <- table$prof_base[liste.na]
      table$prof_sup_moy[liste.na] <- table$prof_sommet[liste.na]
    }
  return(table)
}

lim.hz.NA.2 <- function(table,nom.table){
  # N.B. : fonction adaptée pour repérer les erreurs pour les horizons prof_base et prof_sommet
  if (dim(subset(table, is.na(prof_sommet) | is.na(prof_base)))[1]==0){ 
    print("Pas de NA dans les profondeurs de prélévement")} else {
      
      # 6. cas où prof_base = NA et prof_inf_moy = NA
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- is.na(table$prof_base) & is.na(table$prof_inf_moy)
      mat.na <- subset(table, is.na(table$prof_base) & is.na(table$prof_inf_moy))
      if (nrow(mat.na) !=0){ 
        # identification des derniers horizons 
        liste.end <- vector()
        liste.mid <- vector()
        for (i in 1:nrow(mat.na))
        {
          id.profil.i <-mat.na$id_profil[i]
          hz.i <- mat.na$no_horizon[i]
          liste.hz.i <- subset(table,id_profil == id.profil.i)$no_horizon
          if (hz.i == max(liste.hz.i)) {liste.end <- c(liste.end, mat.na$id_resultat[i])} else {
          liste.mid <- c(liste.mid, mat.na$id_resultat[i])}
        }
        
        mat.na.end <- subset(table, id_resultat %in% liste.end)
        mat.na.mid <- subset(table, id_resultat %in% liste.mid)
        write.table(mat.na.end,paste("lim_hz_NA_2a_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        write.table(mat.na.mid,paste("lim_hz_NA_2b_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        
        table <- subset(table, id_profil %nin% unique(mat.na.mid$id_profil))
        mat.na.end <- subset(mat.na.end , id_profil %nin% unique(mat.na.mid$id_profil))
        table$prof_base[table$id_resultat %in% liste.end] <- table$prof_sommet[table$id_resultat %in% liste.end] + 30
        table$prof_inf_moy[table$id_resultat %in% liste.end] <- table$prof_sommet[table$id_resultat %in% liste.end] + 30
      }
      
    }
  return(table) 
}

### for texture, because we did not extract id_resultat for granulo, I need to use another field, prelevement (I think we calculate the average of the analyses by prelevement....)
lim.hz.NA.2.2 <- function(table,nom.table){
    # N.B. : fonction adaptée pour repérer les erreurs pour les horizons prof_base et prof_sommet
    if (dim(subset(table, is.na(prof_sommet) | is.na(prof_base)))[1]==0){ 
        print("Pas de NA dans les profondeurs de prélévement")} else {
            
            # 6. cas où prof_base = NA et prof_inf_moy = NA
            prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
            colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
            liste.na <- is.na(table$prof_base) & is.na(table$prof_inf_moy)
            mat.na <- subset(table, is.na(table$prof_base) & is.na(table$prof_inf_moy))
            if (nrow(mat.na) !=0){ 
                # identification des derniers horizons 
                liste.end <- vector()
                liste.mid <- vector()
                for (i in 1:nrow(mat.na))
                {
                    id.profil.i <-mat.na$id_profil[i]
                    hz.i <- mat.na$no_horizon[i]
                    liste.hz.i <- subset(table,id_profil == id.profil.i)$no_horizon
                    if (hz.i == max(liste.hz.i)) {liste.end <- c(liste.end, mat.na$id_prelevement[i])} else {
                        liste.mid <- c(liste.mid, mat.na$id_prelevement[i])}
                }
                
                mat.na.end <- subset(table, id_prelevement %in% liste.end)
                mat.na.mid <- subset(table, id_prelevement %in% liste.mid)
                write.table(mat.na.end,paste("lim_hz_NA_2a_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
                write.table(mat.na.mid,paste("lim_hz_NA_2b_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
                
                table <- subset(table, id_profil %nin% unique(mat.na.mid$id_profil))   ### exclude profiles which any horizon is in the middle, and has NA in the base
                mat.na.end <- subset(mat.na.end , id_profil %nin% unique(mat.na.mid$id_profil))  ### exclude profiles which any horizon is in the middle, and has NA in the base
                table$prof_base[table$id_prelevement %in% liste.end] <- table$prof_sommet[table$id_prelevement %in% liste.end] + 30
                table$prof_inf_moy[table$id_prelevement %in% liste.end] <- table$prof_sommet[table$id_prelevement %in% liste.end] + 30
                ### We assume the lower limit of the horizon and prelevement to be 30cm deeper than the upper limit of the prelevement.
            }
            
        }
    return(table) 
}

lim.hz.NA.3 <- function(table,nom.table){
  # N.B. : fonction adaptée pour repérer les horizons où prof_inf_moy = NA mais prof_base != NA
  if (dim(subset(table, is.na(prof_inf_moy) | !is.na(prof_base)))[1]==0){ 
    print("Pas de NA dans les profondeurs de prélévement")} else {
      prof <- as.matrix(cbind(table$prof_sup_moy,table$prof_inf_moy,table$prof_sommet,table$prof_base))
      colnames(prof) <- c("prof_sup_moy","prof_inf_moy","prof_sommet","prof_base")
      liste.na <- !is.na(prof[,"prof_base"]) & is.na(prof[,"prof_inf_moy"]) & (prof[,"prof_sommet"] >= prof[,"prof_sup_moy"])
      mat.na <- subset(table, liste.na)
      if (nrow(mat.na) !=0){ 
        write.table(mat.na,paste("lim_hz_NA_3_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
        table$prof_inf_moy[liste.na] <- table$prof_base[liste.na]
      }
    }
  return(table) 
}

del.na <- function(table, champ, nom.table){
  if (any(is.na(table[,champ])))
  {
    mat.na <- table[is.na(table[,champ]),]
    write.table(mat.na,paste("del_na_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
    table <- subset(table, id_profil %nin% unique(mat.na$id_profil))
  }
  return(table)
}

cor.999 <- function(table, champs, nom.table){
  for (i in 1:length(champs))
  {
    champs.i <- champs[i]
    if (any(table[,champs.i] == 999, na.rm = TRUE))
    {
      mat.999 <- table[table[,champs.i]==999 & !is.na(table[,champs.i]),]
      write.table(mat.999,paste("cor_999_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)    
      table[table[,champs.i]==999 & !is.na(table[,champs.i]),champs.i] <- NA
    }
  }
return(table)
}

# N.B. : il est préférable d'utiliser la fonction cor.topsoil2
cor.topsoil <- function(table, login, password, nom.table){
  ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)
  perso <- odbcConnect(dsn="baracoaPerso",uid=login,pwd=password)
  sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
  
  occup <- sqlQuery(ds3,"
                    SELECT id_profil, occup_codee, occup_libre
                    FROM data.profil")
  
  liste.occup.till <- c(seq(100,508,1),seq(700,1050,1),seq(1200,1400,1), 1600, 3000,5222, seq(6000,6080,1))
  liste.occup.notill <- c(seq(600,603,1),seq(1100,1164,1), 1500, seq(2000,2960,1), seq(4000,5221,1), seq(7000,13100,1))
  
  occup <- subset(occup, id_profil %in% table$id_profil)
  
  # liste des horizons avec  ou non des noms d'horizon
  nom.hz <- unique(subset(table, !is.na(cpcs_nom) | !is.na(rp_2008_nom) | !is.na(rp_95_nom))$id_profil)
  no.nom.hz <- unique(subset(table, is.na(cpcs_nom) & is.na(rp_2008_nom) & is.na(rp_95_nom) & id_profil %nin% nom.hz)$id_profil)
  
  # exportation des horizons à problèmes
  problemes.occup <- subset(occup, id_profil %in% no.nom.hz & !is.na(occup_libre) & is.na(occup_codee))
  if(nrow(problemes.occup) != 0){
    write.table(problemes.occup, "problemes_occup.csv", sep=";", col.names = T,row.names = F) 
  }
    
  # liste des profils labouré et profondeur de labour
  id.prof.till <- vector()
  prof.till.max <- vector() 
  
  # Gestion des profils avec des noms d'horizon renseignés
  if(length(nom.hz) != 0){
    
    # gestion des profils labourés
    for (i in 1:length(nom.hz))
    {
      prof.i <- subset(table, id_profil == nom.hz[i])
      till <- unique(c(grep("L",prof.i$rp_95_nom),grep("Ap",prof.i$rp_95_nom), grep("L",prof.i$rp_2008_nom),grep("p",prof.i$cpcs_nom),grep("P",prof.i$cpcs_nom),grep("LA",prof.i$cpcs_nom)))
      till2 <- unique(c(grep("alu",prof.i$rp_95_nom), grep("alu",prof.i$rp_2008_nom),grep("sal",prof.i$rp_95_nom), grep("sal",prof.i$rp_2008_nom)))
      
      if (length(till)>0 & length(till2)>0){  
        till <- till[till != till2]
      }
      if (length(till)>0){  
        id.prof.till <- c(id.prof.till, unique(prof.i$id_profil))
        hz.till.i <- prof.i[till,]
        prof.till.max <-c(prof.till.max, max(hz.till.i$prof_base))
        prof.sommet.i <- prof.i$prof_sommet
        if (min(prof.sommet.i) !=0){
          id.tot1.i <- hz.till.i$id_tot1[hz.till.i$prof_sommet == min(prof.sommet.i)]
          table$prof_sommet[table$id_tot1 == id.tot1.i] <- 0
        }
      }
    }
    # liste des profils sans labour (selon nom deshorizons)
     id.prof.notill <- setdiff(nom.hz,id.prof.till)
     id.prof.notill <- unique(id.prof.notill)
     no.nom.hz <- c(no.nom.hz, id.prof.notill)
     #if(length(prof.notill) != 0){
     #  prof.i <- subset(table, id_profil == nom.hz[i])
     # if(min(prof.i$prof_sommet) != 0){
     #    table <- subset(table, id_profil != nom.hz[i])
     # }
     #}
  
  }
  
  # Gestion des profils sans nom d'horizon renseigné ou sans nom d'horizon (mais avec une occupation du sol supposant un labour)
  if(length(no.nom.hz) != 0){
 
  occup2 <- subset(occup, id_profil %in% no.nom.hz & !is.na(occup_codee))
  no.nom.hz2 <- no.nom.hz[no.nom.hz %in% occup2$id_profil] # horizons avec occupation du sol codée (à conserver)
  del.nom.hz <- setdiff(no.nom.hz,no.nom.hz2) # horizons avec occupation du sol non codée (à supprimer)
  
  # suppresion des profils sans nom d'horizon, sans occupation du sol et sans horizons de surface
  # on laisse identiques les profils sans nom d'horizon, sans occupation du sol mais avec horizon de surface
  if(length(del.nom.hz) !=0){
      for (i in 1:length(del.nom.hz)){
      prof.i <- subset(table, id_profil == del.nom.hz[i]) 
      if (min(prof.i$prof_sommet) != 0){
        table <- subset(table, id_profil != del.nom.hz[i])
      }
    }
  }
  
  # profils sans nom d'horizon, mais avec occupation du sol
  if(length(no.nom.hz2) != 0){
    for (i in 1:length(no.nom.hz2)){
      # gestion des profils labourés
      if(subset(occup2, id_profil == no.nom.hz2[i])$occup_codee  %in% liste.occup.till){
        prof.i <- orderBy(~no_horizon, subset(table, id_profil == no.nom.hz2[i]))
        prof.sommet.i <- prof.i $prof_sommet
        id.resultat.i <- prof.i$id_resultat[prof.i$prof_sommet == min(prof.sommet.i)]
        # S'il existe un horizon de surface (limite supérieure minimum = 0), on ne fait rien. Sinon : correction.
        # Profils labourés sans horizon de surface, sans horizon décrit dans l'épaisseur labourée. Supression du profil.
        if(min(prof.i$prof_sommet) != 0 & min(prof.i$prof_base) > 35){
          table <- subset(table, id_profil != unique(prof.i$id_profil))
        }
        # Profils labourés sans horizon de surface, avec un horizon décrit dans l'épaisseur labourée. On remonte la limite supérieure minimale à 0.
        if(min(prof.i$prof_sommet)  != 0 & min(prof.i$prof_base) <= 35){
        table$prof_sommet[table$id_resultat == id.resultat.i] <- 0
        }
        if (min(prof.i$prof_sommet)  == 0){
        id.prof.till <- c(id.prof.till, unique(prof.i$id_profil))
        liste.prof.base <- prof.i$prof_base
        id.prof.till.max <- which(liste.prof.base >=25 & liste.prof.base <=35)
        id.prof.till.max <- max(id.prof.till.max)
        if (length(id.prof.till.max)==1){
        prof.till.max <- c(prof.till.max, prof.i$prof_base[id.prof.till.max])  
        } else{
        prof.till.max <- c(prof.till.max, 30)
        }
      }    
    }
  }
 }
}
 till <- as.data.frame(cbind(id.prof.till,rep(1, length(id.prof.till)), prof.till.max))
 till <- rename.vars(till, "V2","tillage")
 table <- merge(table, till, by.x = "id_profil", by.y = "id.prof.till", all.x = T, all.y = F) 
  
odbcCloseAll()
return(table)
}

cor.topsoil2 <- function(table, login, password){
  
  # liste des occupations du sol
  ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)
  # perso <- odbcConnect(dsn="baracoaPerso",uid=login,pwd=password) 
  sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
  
  occup <- sqlQuery(ds3,"
                    SELECT id_profil, occup_codee, occup_libre
                    FROM data.profil")
  
  liste.occup.till <- c(seq(100,508,1),seq(700,1050,1),seq(1200,1400,1), 1600, 3000,5222, seq(6000,6080,1))
  liste.occup.notill <- c(seq(600,603,1),seq(1100,1164,1), 1500, seq(2000,2960,1), seq(4000,5221,1), seq(7000,13100,1))
  
  occup <- subset(occup, id_profil %in% table$id_profil)
  
  # liste des profils
  liste.profil <- unique(table$id_profil)
  
  # barre indiquant la durée du processus
  pb <- txtProgressBar(min = 0, max = length(liste.profil) , style = 3)
  
  # initialisation des listes de profils labourés et des profondeurs de labour
  id.prof.till <- vector()
  prof.till.max <- vector()
  
  for (i in 1:length(liste.profil)){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    
    table.i <- subset(table, id_profil == liste.profil[i])
    id.prof.till.i <- vector()
    prof.till.max.i <- vector()
    
    # profil avec des nom d'horizon renseignés pour tous les horizons du profil
    if (FALSE %nin% unique(!is.na(table.i$cpcs_nom) | !is.na(table.i$rp_2008_nom) | !is.na(table.i$rp_95_nom))){
      
      table.i$rp_95_nom <- tolower(table.i$rp_95_nom)
      table.i$rp_2008_nom <- tolower(table.i$rp_2008_nom)
      table.i$cpcs_nom <- tolower(table.i$cpcs_nom)
      
      till <- unique(c(grep("l",table.i$rp_95_nom),grep("ap",table.i$rp_95_nom), grep("l",table.i$rp_2008_nom),grep("p",table.i$cpcs_nom),grep("la",table.i$cpcs_nom)))
      till2 <- unique(c(grep("bp",table.i$cpcs_nom),grep("jp",table.i$cpcs_nom),grep("alu",table.i$rp_95_nom), grep("alu",table.i$rp_2008_nom),grep("sal",table.i$rp_95_nom), grep("sal",table.i$rp_2008_nom)))
            
      if (length(till)>0 & length(till2)>0){  
        till <- setdiff(till, till2)
      }
      if (length(till)>0){  
        id.prof.till.i <- unique(table.i$id_profil)
        hz.till.i <- table.i[till,]
        hz.till.i <- subset(hz.till.i, prof_base >= 25 & prof_base <=35)
        if (nrow(hz.till.i) >= 1){
          prof.till.max.i <- max(hz.till.i$prof_base)
        } else {
          prof.till.max.i <- 30 
        }
        prof.sommet.i <- table.i$prof_sommet
        if (min(prof.sommet.i) !=0){
          id.tot1.i <- hz.till.i$id_tot1[hz.till.i$prof_sommet == min(prof.sommet.i)]
          table.i$prof_sommet[table.i$id_tot1 == id.tot1.i] <- 0
        }
      }
      
    }
    
    if (FALSE %in% unique(!is.na(table.i$cpcs_nom) | !is.na(table.i$rp_2008_nom) | !is.na(table.i$rp_95_nom)) | length(prof.till.max.i)==0){
      # profil sans nom d'horizon (total ou partiel)  
      occup.i <- subset(occup, id_profil == unique(table.i$id_profil))$occup_codee
      
      # si il n'y a ni nom d'horizon, ni occupation du sol, et que l'horizon supérieur n'a pas 0 en limite supérieure : suppression du profil
      if (length(occup.i) != 1){
        if (min(table.i$prof_sommet) != 0){
          next
        }        
      }
      
      # si occupation du sol correspond à une parcelle labourée
      if (occup.i  %in% liste.occup.till){
        table.i <- orderBy(~no_horizon, table.i)
        prof.sommet.i <- table.i $prof_sommet
        id.resultat.i <- table.i$id_resultat[table.i$prof_sommet == min(prof.sommet.i)]
        # S'il existe un horizon de surface (limite supérieure minimum = 0), on ne fait rien. Sinon : correction.
        # Profils labourés sans horizon de surface, sans horizon décrit dans l'épaisseur labourée. Supression du profil.
        if(min(table.i$prof_sommet) != 0 & min(table.i$prof_base) > 35){
          next
        }
        # Profils labourés sans horizon de surface, avec un horizon décrit dans l'épaisseur labourée. On remonte la limite supérieure minimale à 0.
        if(min(table.i$prof_sommet)  != 0 & min(table.i$prof_base) <= 35){
          table.i$prof_sommet[table.i$id_resultat == id.resultat.i] <- 0
        }
        if (min(table.i$prof_sommet)  == 0){
          id.prof.till.i <- unique(table.i$id_profil)
          liste.prof.base <- table.i$prof_base
          id.prof.till.max <- which(liste.prof.base >=25 & liste.prof.base <=35)
          if(length(id.prof.till.max) == 1){
            id.prof.till.max <- max(id.prof.till.max)
            prof.till.max.i <- table.i$prof_base[id.prof.till.max] 
          }else{
            prof.till.max.i <- 30
          }
        }  
      } 
    }
    
    prof.till.max <- c(prof.till.max, prof.till.max.i)
    id.prof.till <- c(id.prof.till, id.prof.till.i)
    if (i == 1){
      table.r <- table.i
    } else {
      table.r <- rbind(table.r, table.i)
    }
  }
  
  close(pb)
  
  till <- as.data.frame(cbind(id.prof.till,rep(1, length(id.prof.till)), prof.till.max))
  till <- rename.vars(till, "V2","tillage", info= F)
  table.r <- merge(table.r, till, by.x = "id_profil", by.y = "id.prof.till", all.x = T, all.y = F) 
  
  odbcCloseAll()
  return(table.r)    
}

#### Because I work with the upper and lowr horizon limits (prof_sup_moy and prof_inf_moy), I change these names in the function.
#### for texture, because we did not extract id_resultat for granulo, I need to use another field, prelevement (I think we calculate the average of the analyses by prelevement....)
cor.topsoil3 <- function(table, login, password){
    
    # liste des occupations du sol
    ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)
    # perso <- odbcConnect(dsn="baracoaPerso",uid=login,pwd=password) 
    sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 
    
    occup <- sqlQuery(ds3,"
                      SELECT id_profil, occup_codee, occup_libre
                      FROM data.profil")
    
    liste.occup.till <- c(seq(100,508,1),seq(700,1050,1),seq(1200,1400,1), 1600, 3000,5222, seq(6000,6080,1))
    liste.occup.notill <- c(seq(600,603,1),seq(1100,1164,1), 1500, seq(2000,2960,1), seq(4000,5221,1), seq(7000,13100,1))
    
    occup <- subset(occup, id_profil %in% table$id_profil)
    
    # liste des profils
    liste.profil <- unique(table$id_profil)
    
    # barre indiquant la durée du processus
    pb <- txtProgressBar(min = 0, max = length(liste.profil) , style = 3)
    
    # initialisation des listes de profils labourés et des profondeurs de labour
    id.prof.till <- vector()
    prof.till.max <- vector()
    
    for (i in 1:length(liste.profil)){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, i)
        
        
        table.i <- subset(table, id_profil == liste.profil[i])
        id.prof.till.i <- vector()
        prof.till.max.i <- vector()
        
        # profil avec des nom d'horizon renseignés pour tous les horizons du profil
        if (FALSE %nin% unique(!is.na(table.i$cpcs_nom) | !is.na(table.i$rp_2008_nom) | !is.na(table.i$rp_95_nom))){
            
            table.i$rp_95_nom <- tolower(table.i$rp_95_nom)
            table.i$rp_2008_nom <- tolower(table.i$rp_2008_nom)
            table.i$cpcs_nom <- tolower(table.i$cpcs_nom)
            
            till <- unique(c(grep("l",table.i$rp_95_nom),grep("ap",table.i$rp_95_nom), grep("l",table.i$rp_2008_nom),grep("p",table.i$cpcs_nom),grep("la",table.i$cpcs_nom)))
            till2 <- unique(c(grep("bp",table.i$cpcs_nom),grep("jp",table.i$cpcs_nom),grep("alu",table.i$rp_95_nom), grep("alu",table.i$rp_2008_nom),grep("sal",table.i$rp_95_nom), grep("sal",table.i$rp_2008_nom)))
            
            if (length(till)>0 & length(till2)>0){  
                till <- setdiff(till, till2)
            }
            if (length(till)>0){  
                id.prof.till.i <- unique(table.i$id_profil)
                hz.till.i <- table.i[till,]
                hz.till.i <- subset(hz.till.i, prof_inf_moy >= 25 & prof_inf_moy <=35)
                if (nrow(hz.till.i) >= 1){
                    prof.till.max.i <- max(hz.till.i$prof_inf_moy)
                } else {
                    prof.till.max.i <- 30 
                }
                prof.sup.i <- table.i$prof_sup_moy
                if (min(prof.sup.i) !=0){
                    id.tot1.i <- hz.till.i$id_tot1[hz.till.i$prof_sup_moy == min(prof.sup.i)]
                    table.i$prof_sup_moy[table.i$id_tot1 == id.tot1.i] <- 0
                }
            }
            
        }
        
        if (FALSE %in% unique(!is.na(table.i$cpcs_nom) | !is.na(table.i$rp_2008_nom) | !is.na(table.i$rp_95_nom)) | length(prof.till.max.i)==0){
            # profil sans nom d'horizon (total ou partiel)  
            occup.i <- subset(occup, id_profil == unique(table.i$id_profil))$occup_codee
            
            # si il n'y a ni nom d'horizon, ni occupation du sol, et que l'horizon supérieur n'a pas 0 en limite supérieure : suppression du profil
            if (length(occup.i) != 1){
                if (min(table.i$prof_sup_moy) != 0){
                    next
                }        
            }
            
            # si occupation du sol correspond à une parcelle labourée
            if (occup.i  %in% liste.occup.till){
                table.i <- orderBy(~no_horizon, table.i)
                prof.sup.i <- table.i $prof_sup_moy
                id.prelevement.i <- table.i$id_prelevement[table.i$prof_sup_moy == min(prof.sup.i)]
                # S'il existe un horizon de surface (limite supérieure minimum = 0), on ne fait rien. Sinon : correction.
                # Profils labourés sans horizon de surface, sans horizon décrit dans l'épaisseur labourée. Supression du profil.
                if(min(table.i$prof_sup_moy) != 0 & min(table.i$prof_inf_moy) > 35){
                    next
                }
                # Profils labourés sans horizon de surface, avec un horizon décrit dans l'épaisseur labourée. On remonte la limite supérieure minimale à 0.
                if(min(table.i$prof_sup_moy)  != 0 & min(table.i$prof_inf_moy) <= 35){
                    table.i$prof_sup_moy[table.i$id_prelevement == id.prelevement.i] <- 0
                }
                if (min(table.i$prof_sup_moy)  == 0){
                    id.prof.till.i <- unique(table.i$id_profil)
                    liste.prof.inf <- table.i$prof_inf_moy
                    id.prof.till.max <- which(liste.prof.inf >=25 & liste.prof.inf <=35)
                    if(length(id.prof.till.max) == 1){
                        id.prof.till.max <- max(id.prof.till.max)
                        prof.till.max.i <- table.i$prof_inf_moy[id.prof.till.max] 
                    }else{
                        prof.till.max.i <- 30
                    }
                }  
            } 
        }
        
        prof.till.max <- c(prof.till.max, prof.till.max.i)
        id.prof.till <- c(id.prof.till, id.prof.till.i)
        if (i == 1){
            table.r <- table.i
        } else {
            table.r <- rbind(table.r, table.i)
        }
    }
    
    close(pb)
    
    till <- as.data.frame(cbind(id.prof.till,rep(1, length(id.prof.till)), prof.till.max))
    till <- rename.vars(till, "V2","tillage", info= F)
    table.r <- merge(table.r, till, by.x = "id_profil", by.y = "id.prof.till", all.x = T, all.y = F) 
    
    odbcCloseAll()
    return(table.r)    
}


MO.to.C <- function(table){
  table$valeur[!is.na(table$no_methode) & table$no_methode == 29] <- table$valeur[!is.na(table$no_methode) & table$no_methode == 29] / 2
  table$valeur[!is.na(table$no_methode) & table$no_methode != 29] <- table$valeur[!is.na(table$no_methode) & table$no_methode != 29]/ 1.724
  table$valeur[is.na(table$no_methode)] <- table$valeur[is.na(table$no_methode)]/ 1.724 
  return(table)
}

select.seuil <- function(table, seuil,critere){
  # critère = 0 --> sélection des profils avec horizons strictement inférieurs au seuil
  # critère = 1 --> sélection des profils avec horizons strictement inférieurs au seuil
  # critère = 2 --> sélection des profils avec horizons inférieurs ou égal au seuil 
  # critère = 3 --> sélection des profils avec horizons strictement supérieurs au seuil
  # critère = 4 --> sélection des profils avec horizons supérieurs ou égal au seuil 
  
  if (critere == 0){
    hz <- unique(subset(table, valeur != seuil)$id_profil)
    table <- subset(table, id_profil %nin% hz)
  }
  if (critere == 1){
    hz <- unique(subset(table, valeur >=  seuil)$id_profil)
    table <- subset(table, id_profil %nin% hz)
  }
    if (critere == 2){
      hz <- unique(subset(table, valeur > seuil)$id_profil)
      table <- subset(table, id_profil %nin% hz)
  }
  if (critere == 3){
    hz <- unique(subset(table, valeur <= seuil)$id_profil)
    table <- subset(table, id_profil %nin% hz)
  }
  if (critere == 4){
    hz <- unique(subset(table, valeur < seuil)$id_profil)
    table <- subset(table, id_profil %nin% hz)
  }
  return(table)
}

homogenize <- function(table1, table2, output){
  if ((ncol(table1) != ncol(table2)) ){
    stop("Nombre de colonnes différent")
  }
  if(length(intersect(colnames(carbone), colnames(MO))) != ncol(table1)){
    stop("Noms de colonne différents")
  }
  nom.col <- colnames(table1)
  for (i in 1: length(nom.col)){
    if (is.factor(table1[,nom.col[i]]) | is.factor(table2[,nom.col[i]])) {
      #liste.levels <- union(levels(table1[,nom.col[i]]), levels(table2[,nom.col[i]]))
      #levels(table1[,nom.col[i]]) <-  liste.levels 
      #levels(table2[,nom.col[i]]) <-  liste.levels
      table1[,nom.col[i]] <- as.character(table1[,nom.col[i]])
      table2[,nom.col[i]] <- as.character(table2[,nom.col[i]])
    }
  }
  for (i in 1: length(nom.col)){
    if (class(table1[,nom.col[i]]) != class(table2[,nom.col[i]])) {
      stop(paste( nom.col[i]," : type de colonne différent. Voir fonction class()"))
     }
  }
  if(output == 1) return(table1)
  if(output == 2) return(table2)
}

profsol.cor <- function(table, champ.id, champ.prof.sup, champ.prof.inf, champ.noms.hz, champ.prof.sol){
  # DESCRIPTION : fonction permettant de sélectionner uniquement les horizons de sol (supressions des horizons C, R, M et D)
  #               et de renseigner la variable "profondeur de sol" si possible (à partir du nom des horizons s'ils sont renseignés)
  # INPUTS
  
  liste <- unique(table[,champ.id])
  pb <- txtProgressBar(min = 0, max = length(liste) , style = 3)
  
  for (i in 1:length(liste)){
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    
    prof.i <- table[table[,champ.id]== liste[i],]
    psol.i <- unique(prof.i[,champ.prof.sol])
    
    if (is.na(psol.i)){     # cas où la profondeur du sol n'est pas renseignée
      
      liste.hz.i <- rep(NA, nrow(prof.i))
      
      for (j in 1:length(champ.noms.hz)){
        liste.hz.i2 <- as.character(prof.i[,champ.noms.hz[j]]) 
        liste.hz.i[is.na(liste.hz.i)] <- liste.hz.i2[is.na(liste.hz.i)] 
      }
      
      pmat.i <- unique(c(grep("C", liste.hz.i),grep("R", liste.hz.i), grep("M",prof.i$liste.hz.i)))
      
      # si un horizon C, R, M ou D est décrit, on donne la limite supérieure de l'horizon comme profondeur du sol
      # sinon on ne fait rien
      if (length(pmat.i)>0){
        pmat.i <- min(pmat.i)
        psol.i2 <- prof.i[pmat.i,champ.prof.sup]
        prof.i[,champ.prof.sol] <- psol.i2 
        prof.i <- subset(prof.i, prof_sup_moy < psol.i2)
        
      }
      
    } else {                 # cas où la profondeur du sol est renseignée
      
      prof.i <- subset(prof.i, prof_sup_moy < psol.i)
    }
    
    if (i == 1){
      prof.r <- prof.i
    } else {
      prof.r <- rbind(prof.r,prof.i)
    }
    
  }
  
  close(pb)
  return(prof.r)  
}

cor.string <- function(x){
  x <- gsub("’","",x)
  x <- gsub("'","",x)
  x <- gsub("\"","",x)
  x <- gsub("\\.","",x)
  x <- gsub("Ã©","e",x)
  x <- gsub("Ã¨","e",x)
  x <- gsub("Ãš¨","e",x)
  x <- gsub("Ã¨","a",x)
  x <- gsub("Ã“","a",x)
  x <- gsub("Ã","a",x)
  x <- gsub("Ã´","o",x)
  x <- gsub("  "," ",x)
  x <- gsub("  "," ",x)
  x <- gsub("  "," ",x)
  x <- tolower(x)
  x <- gsub("ã","a",x)
  x <- gsub("aš","e",x)
  x <- gsub("é","e",x)
  x <- gsub("è","e",x)
  x <- gsub("ê","e",x)
  x <- gsub("ë","e",x)
  x <- gsub("à","a",x)
  x <- gsub("â","a",x)
  x <- gsub("ô","o",x)
  x <- gsub("ï","i",x)
  x <- gsub("î","i",x)
  x <- gsub("ù","u",x)
  x <- gsub("û","u",x)
  x <- gsub("ü","u",x)
  x <- gsub("œ","oe",x)
  x <- gsub("a‰","e",x)
  x <- gsub("aˆ","e",x)
  x <- gsub("až","i",x)
  x <- gsub("a€","a",x)
  x <- gsub("a¯","i",x)
  x <- gsub("gr+s","gres",x)
  x <- gsub("^\\s+|\\s+$", "",x)
  return(x)
}

replace <- function(data){
  for(i in 1:ncol(data)){
    
    if(is.factor(data[,i])){
      m <- as.matrix(data[,i])
      for(j in 1:nrow(data)){
        m[j,1] <- cor.string(m[j,1])
      }
      data[,i] <- m
    }
  }
  return(data)
} 

select.mot1 <- function(x){
  x <- unlist(strsplit(x, " "))[1]
}

add.id <- function(table, champs, nom.id){
  new.id <- data.matrix(table[,champs],rownames.force=FALSE)
  new.id <- apply(new.id, 1, conc.champ)
  table$new.id  <- new.id 
  table <- rename.vars(table, "new.id",nom.id, info=FALSE)
  return(table)
}

MPE <- function(x, y){
  return (mean(x-y))
}  

SDPE <- function(x, y){
  return ((sd(x - y)))
}

RMSPE <- function(x, y){
  return (mean((x-y)^2)^0.5)
}

MedPE <- function(x, y){
  return (median(x-y))
}	

RMedSPE <- function(x, y){
  return (median((x-y)^2)^0.5)
}

R2 <- function(x, y){
  return (cor(x,y)^2)
}

summaryStats <- function(pred, obs){
  #c(paste( "MPE : ", MPE(pred, obs)),
  #   paste( "SDPE : ", SDPE(pred, obs)),
  #   paste( "RMSPE : ", RMSPE(pred, obs)),
  #   cor(pred, obs)^2)
  res <- c(cor(pred, obs)^2, MPE(pred, obs), SDPE(pred, obs), RMSPE(pred, obs))
  names(res) <- c("R2", "MPE", "SDPE", "RMSPE")
  return (res)
  #   
}

jeffrey <- function(X, a1, a2) {
  C <- X	
  expr <- expression(a1 + a2 * log10(C))
  eval(expr)
}

federer <- function(X, a1, a2, a3) {
  C <- X	
  expr <- expression(a1 + a2 * log(C) + a3 * log(C)^2)
  eval(expr)
}

adams <- function(X, a1, a2) {
  C <- X	
  expr <- expression(100/(a1 + a2 * C))
  eval(expr)
}

manrique <- function(X, a1, a2) {
  C <- X	
  expr <- expression(a1 + a2 * C^0.5)
  eval(expr)
}

kaur <- function(X, a1, a2, a3, a4, a5) {
  C <- X[,1]
  a <- X[,2]
  lim <- X[,3] + X[,4]
  expr <- expression(a1 + a2 * C + a3 * a + a4 * a^2 + a5 * lim)
  eval(expr)
}

