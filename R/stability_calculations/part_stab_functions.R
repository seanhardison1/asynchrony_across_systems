# This script was adapted for our work and the original can be found in:
# 
# Data S1 to Lamy T, Wang S, Renard D, Lafferty KD, Reed DC and Miller RJ.
# Species insurance trumps spatial insurance in stabilizing biomass of a marine macroalgal metacommunity. Ecology.
# 
# As part of:
# 
# Lamy, T., Wang, S., Renard, D., Lafferty, K. D., Reed, D. C., & Miller, R. J. (2019). 
# Species insurance trumps spatial insurance in stabilizing biomass of a marine macroalgal metacommunity. 
# Ecology, 100(7), e02719.
# 
# ------
# 
# File contains two R functions, part_stab_comp() and part_stab_test(),
# to investigate ecological variability across hierarchical levels
#


part_stab_comp <- function(Y, s, t,
                           Y2 = NULL, s2 = NULL, t2 = NULL,
                           surv_weight = F, cv = F)
  #
  # Description --
  # R function to partition ecological variability across hierarchical levels
  # The function partitions ecological variability across two spatial scales (local and regional) 
  # and two organizational levels (population and community)
  # Six synchrony metrics serve as scaling factors to measure how variability scales from one specific hierarchical level to the next 
  # These synchrony metrics are the square-root transformation of the synchrony metric developed by Loreau and de Mazancourt (2008)
  #
  # Arguments --
  #
  # Y: community table: an observation by species matrix with row blocks corresponding to local patches
#    contains the biomass of species (i) in local patch (k) 
# s: number of patches
# t: number of years
#
# Outputs --
# Overall table:
# Species-level variability at the local scale = CV_SL
# Community-level variability at the local scale = CV_CL
# Species-level variability at the regional scale = CV_SR (weighted average of metapopulation variability across species)
# Community-level variability at the regional scale = CV_CL (variability of the whole metacommunity)
# Average species-level spatial synchrony = phi_S_LR
# Community-level spatial synchrony = phi_C_LR
# Average local-scale species synchrony = phi_SC_L
# Regional-scale species synchrony = phi_SC_R
#
# species-level table (i):
# Spatial synchrony of species i = phi_i_LR
# weight of species i = w_i
#
# pacth-level table (k):
# Species synchrony within patch k = phi_SC_k
# weight of pacth k = w_k
#
# Reference --
#
# Wang S, Lamy T, Hallett LM and Loreau M (2019) Stability and synchrony across ecological hierarchies in heterogeneous metacommunities: linking theory to data.
# Ecography. doi: 10.1111/ecog.04290.
# Wang and Loreau (2014) Ecosystem stability in space: alpha, beta and gamma variability.
# Ecology Letters. 17: 891-901.
# Loreau M and de Mazancourt C (2008) Species synchrony and its drivers: neutral and nonneutral community dynamics in fluctuating environments.
# American Naturalist. 172: E48-66.
# 
# Author: Thomas LAMY
# edited September 28th 2016
{
  # Total community biomass (sum within rows)
  tot.bio <- apply(Y, 1, sum, na.rm = T)
 
  
  # Check that data is balanced
  if(length(tot.bio) != s*t) stop("STOP: sites are not surveyed every year")
  
  
  # Matrix of total community biomass within each patch
  Y.CL <- matrix(tot.bio, nrow=t, ncol=s, byrow=FALSE)
  
  
  # List of local patches containing population biomass of each species
  # for each of 20 sites
  SiteL <- list()
  for(k in 1:s) {
    SiteL[[k]] <- Y[c(((k-1)*t+1):(k*t)),]
  }

  
  # Matrix of metapopulation biomass
  Y.SR <- Reduce("+", SiteL) 
  
  
  # Temporal mean biomass of the whole metacommunity
  mu_TT <- mean(apply(Y.SR, 1, sum, na.rm = T))
  
  
  # Temporal sd of species i in patch k
  sd_ik <- matrix(unlist(lapply(SiteL, function(x){apply(x, 2, sd, na.rm = T)})), 
                  nrow=s, ncol=dim(Y)[2], byrow=TRUE)

  
  # Temporal mean biomass of species i in patch k
  mu_ik <- matrix(unlist(lapply(SiteL, function(x){apply(x, 2, mean, na.rm = T)})), 
                  nrow=s, ncol=dim(Y)[2], byrow=TRUE)

  
  # Temporal cv of species i in patch k
  cv_ik <- sd_ik/mu_ik
  
  
  # Temporal mean total community biomass in patch k
  mu_Tk <- matrix(unlist(lapply(SiteL, function(x){mean(apply(x, 1, sum, na.rm = T))})), 
                  nrow=1, ncol=s, byrow=TRUE)

  
  # Temporal sd of total community biomass in patch k
  sd_Tk <- matrix(unlist(lapply(SiteL, function(x){sd(apply(x, 1, sum, na.rm = T))})), 
                  nrow=1, ncol=s, byrow=TRUE)

 
  # Temporal cv of total community biomass in patch k
  cv_Tk <- (sd_Tk/mu_Tk) 

  
  # temporal sd of species i at the regional scale
  sd_iT <- apply(Y.SR, 2, sd, na.rm = T)

  # temporal mean biomass of species i at the regional scale
  mu_iT <- apply(Y.SR, 2, mean, na.rm = T)

  # temporal cv of species i at the regional scale
  cv_iT <- sd_iT/mu_iT

  # Species-level variability within patch k
  # defined as the weighted average of species variability across species
  cv_SL_k <- rowSums(sd_ik)/mu_Tk

  
  # Average species-level variability at the local scale
  # defined as the weighted average of species-level variability across patches
  CV_SL_part <- data.frame(CV_SL_part = as.numeric(cv_SL_k),
                           CV_SL_weight = as.numeric(mu_Tk/mu_TT))
  CV_SL <- sum(mu_Tk/mu_TT * cv_SL_k)
  
  # Average community-level variability at the local scale
  # defined as the weighted average of community variability across patches
  CV_CL_part <- data.frame(CV_CL_part = as.numeric(cv_Tk),
                           CV_CL_weight = as.numeric(mu_Tk/mu_TT),
                           sd_Tk = as.numeric(sd_Tk),
                           mu_Tk = as.numeric(mu_Tk))
  CV_CL <- sum(mu_Tk/mu_TT * cv_Tk)
  
  
  # print(cv_Tk)
  # Species-level variability at the regional scale
  # defined as the weighted average of metapopulation variability across species
  CV_SR_part <- data.frame(CV_SR_part = cv_iT,
                           SD_SR_part = sd_iT,
                           weight = mu_iT/mu_TT)
  CV_SR <- sum(mu_iT/mu_TT * cv_iT)
  SD_SR <- sum(mu_iT/mu_TT * sd_iT)

  # Community-level variability at the regional scale, or the variability of the whole metacommunity
  # variance-covariance matrix of the temporal biomass across sites
  W <- cov(Y.CL)
  
  # Sum of the temporal covariances of the communities in each patch
  Covsum <- sum(W)
  
  # Temporal variability at the metacommunity scale is the coefficient of temporal variation of metacommunity biomass
  CV_CR <- sqrt(Covsum)/mu_TT
  SD_CR <- sqrt(Covsum)
  
  # Species synchrony within patch k
  phi_SC_k <- as.numeric(sd_Tk/rowSums(sd_ik))
  
  # Average local-scale species synchrony
  # defined as the weighted average of species synchrony across patches
  w_k <- rowSums(sd_ik)/sum(sd_ik)
  
  phi_SC_L <- sum(w_k * phi_SC_k)
  
  # Spatial synchrony of species i
  phi_i_LR <- sd_iT/colSums(sd_ik)
  
  # Average species-level spatial synchrony
  # defined as the weighted average of spatial synchrony across species
  w_i <- colSums(sd_ik)/sum(sd_ik) 
  phi_S_LR <- sum(w_i * phi_i_LR, na.rm=TRUE)
  
  # Community-level spatial synchrony
  phi_C_LR <- sqrt(Covsum)/((sum(sd_Tk)))
  
  # Regional-scale species synchrony
  phi_SC_R <- sqrt(Covsum)/sum(sd_iT)
  
  # Overall table
  part <- data.frame(CV_SL = CV_SL, 
                     CV_CL = CV_CL, 
                     CV_SR = CV_SR,
                     SD_SR = SD_SR,
                     SD_CR = sqrt(Covsum),
                     CV_CR = CV_CR,
                     SD_CR = SD_CR,
                     mu_TT = mu_TT,
                     phi_SC_L=phi_SC_L, phi_S_LR=phi_S_LR, 
                     phi_C_LR=phi_C_LR, phi_SC_R=phi_SC_R)
  # Species-level table (i)
  spe <- data.frame(phi_i_LR=phi_i_LR, w_i=w_i)
  # Patch-level table (k)
  if (cv){
    patch <- data.frame(phi_SC_k=phi_SC_k, w_k=w_k,
                        cv_Tk = as.vector(cv_Tk))
  } else {
    patch <- data.frame(phi_SC_k=phi_SC_k, w_k=w_k)
  }
  
  res <- list(part=part, spe=spe, patch=patch, CV_SR_part = CV_SR_part, 
              CV_CL_part = CV_CL_part,
              CV_SL_part = CV_SL_part)
  res
}

part_stab_test <- function(Y, s, t, nrands=999, fit_gams = T,
                           fit_trend = T)
  #
  # Description --
  # R function to test the synchrony metrics across hierarchical levels
  # Each synchrony metric serves as scaling factors to measure how variability scales across two spatial scales (local and regional) 
  # and two organizational levels (population and community) 
  # These synchrony metrics are the square-root transformation of the synchrony metric developed by Loreau and de Mazancourt (2008):
  # phi=sigma_X/(sigma_x), whith sigma_X the temporal standard deviation of biomass at the higher hierarchical level X, 
  # and sigma_x the sum of the temporal standard deviation of biomass at the lower hierarchical levels x.
  #
  # The significance of all synchrony indices is tested by generating null distributions based on cyclic-shift permutations
  # For each synchrony index, the permutations maintain the observed temporal variability and temporal autocorrelation of the lower hierarchical levels x
# but generate random distributions for the higher hierarchical level X
# given either spatially independent fluctuations across:
# local populations of a given species (phi_i_LR and phi_S_LR)
# local communities (phi_C_LR)
# or given independent fluctuations among:
# local populations within a given patch (phi_SC_k and phi_SC_L)
# metapopulations (phi_SC_R) 
#
# Arguments --
#
# Y: community table: an observation by species matrix with row blocks corresponding to local patches
#    contains the biomass of species (i) in local patch (k) 
# s: number of patches
# t: number of years
# nrands: the number of randomizations
#
# Outputs --
# Overall table:
# Average species-level spatial synchrony = phi_S_LR
# Community-level spatial synchrony = phi_C_LR
# Average local-scale species synchrony = phi_SC_L
# Regional-scale species synchrony = phi_SC_R
#
# species-level table (i):
# Spatial synchrony of species i = phi_i_LR
#
# pacth-level table (k):
# Species synchrony within patch k = phi_SC_k
#
# For each of these three tables this function returns the following metrics:
# pval_greater and pval_lower correspond to the probability that a given synchrony metric is higher or lower, respectively, than expected given the null model
# lw95 and up95 correspond to the 95% confidence interval of the randomizations
# SES corresponds to the Standardized Effect Size obtained from the null distribution
# SES = (phi-mean_null)/sd_null, with mean_null and sd_null the mean and standard deviation of the null distribution
# perc is the percentage by which a given synchrony metric is higher or lower than expected given the randomizations
#
# Reference --
#
# Wang S, Lamy T, Hallett LM and Loreau M (2019) Stability and synchrony across ecological hierarchies in heterogeneous metacommunities: linking theory to data.
# Ecography. doi: 10.1111/ecog.04290.
# Loreau M and de Mazancourt C (2008) Species synchrony and its drivers: neutral and nonneutral community dynamics in fluctuating environments.
# American Naturalist. 172: E48-66.
# Hallett LM, et al. (2014). Biotic mechanisms of community stability shift along a precipitation gradient. Ecology. 95: 1693-1700.
# 
# Author: Thomas LAMY
# edited December 18th 2018
{
  #### Internal functions
  ## General function to compute synchrony metrics following the definition of
  ## Loreau & de Mazancourt (2008) but in a square root version
  sync.fun <- function(Y){
    # temporal variance of hierarchical level X 
    X.var <- var(rowSums(Y))
    # standard deviations of the lower hierarchical levels (x) making up the hierarchical level X
    x.sd <- apply(Y, 2, sd)
    return(sqrt(X.var/sum(x.sd, na.rm=TRUE)^2))
  }
  
  ## function to perform cyclic-shift permutation as adapted from Hallett et al. (2014)
  cs_perm <- function(Y, ng = NULL){
    
    Y <- as.matrix(Y)
    
    nr <- nrow(Y);  nc <- ncol(Y)
    rand.Y <- matrix(NA, nrow=nr, ncol=nc)
    for (i in seq_len(nc)) rand.Y[,i] <- permute::shuffleSeries(Y[,i])
    # print(matplot(rand.Y))
    rand.Y <- as.data.frame(rand.Y)
    return(rand.Y)
  }
  
  # cs_perm(1:10)
  
  # cs_perm2(Y = sflo2$est)
  
  mod_fun <- function(df, fit_trend = T){
    tryCatch(
      expr = {
      # fit trend model and extract coefficient
      if (fit_trend){
        mod <- gam(value ~ s(year),
                   family = "tw", data = df) 
      } else if (!fit_trend) {
        mod <- gam(value ~ 1,
                   family = "tw", data = df)
      }
      
      mu <- predict(mod, type = "response")
      p <- as.numeric(str_extract(mod$family[[1]][1], "1\\.\\d{1,3}"))
      
      # model residual correlation structure
      residual_arima <- 
        auto.arima(residuals.gam(mod),
                   stationary = TRUE, seasonal = FALSE)
      residual_arima$p <- p
      residual_arima$mu <- mu
      return(residual_arima)
        },
      error = function(e){
        print("Model failed.")
      return(NA)
      }
    )
  }
  
  cs_perm2 <- function(Y, ng, series, fit_trend = T){
    mod_set <- 
      Y %>%
      as.data.frame() %>% 
      gather(var, value) %>% 
      mutate(year = rep(1:series, ng)) %>% 
      group_by(var) %>% 
      nest() %>% 
      mutate(sim = map(data,
                       mod_fun,
                       fit_trend = fit_trend))
    return(mod_set)
  }
  
  # Total community biomass
  tot.bio <- apply(Y, 1, sum)
  # check that data is balanced
  if(length(tot.bio) != s * t) stop("STOP: sites are not surveyed every years")
  # Matrix of total community biomass within each patch
  Y.CL <- matrix(tot.bio, nrow=t, ncol=s, byrow=FALSE)
  # List of local communities containing population biomass of each species
  SiteL <- list(); for(k in 1:s) SiteL[[k]] <- Y[c(((k-1)*t+1):(k*t)),] 
  # Matrix of metapopulation biomass
  Y.SR <- Reduce("+", SiteL) 
  
  # Partitioning ecological stability across hierarchical levels using part_stab_comp
  stab <- part_stab_comp(Y, s, t)
  obs <- stab$part
  obs_spe <- as.matrix(stab$spe)
  obs_patch <- as.matrix(stab$patch)
  
  # Testing regional-scale species synchrony (phi_S->C,R)  
  phi_SC_R_rands <- numeric(length=nrands+1)*NA
  cat("Testing regional-scale species synchrony \n")
  prog.bar=txtProgressBar(min=0, max=nrands, style=3)
  if (fit_gams){
    mods <- cs_perm2(Y.SR, ng = ncol(Y.SR), series = nrow(Y.SR),
                     fit_trend)
    sim_out <- NULL
    for (j in 1:nrow(mods)){
      # extract time series model for each species
      m <- 
        mods %>% 
        ungroup() %>% 
        slice(j) %>% 
        .$sim
      
      if (!all(is.na(m[[1]]))){
        p <- m[[1]]$p
        mu <- m[[1]]$mu
        
        sim <-  replicate(n = nrands,
                          rTweedie(mu, p) + 
                            exp(as.numeric(simulate(m[[1]])))
        ) %>% 
          as.data.frame() %>% 
          gather(var, value) %>% 
          mutate(grp = j)
        
      } else {
        sim <- tibble(var = rep(paste0("V",1:nrands), each = nrow(mods$data[[1]])),
                      value = 0,
                      grp = j)
      }
      
      
      assign("sim_out", rbind(sim_out, sim))
      setTxtProgressBar(prog.bar, j)
    }
    
    # simulated_metapops <-
    #   ggplot(sim_out %>%
    #            group_by(grp, var) %>%
    #            mutate(year = 2002:2018,
    #                   common = case_when(grp == 1 ~ "Atlantic croaker",
    #                                      grp == 2 ~ "black drum",
    #                                      grp == 3 ~ "spot",
    #                                      grp == 4 ~ "striped bass",
    #                                      grp == 5 ~ "summer flounder",
    #                                      grp == 6 ~ "weakfish",
    #                                      grp == 7 ~ "white perch"))) +
    #   geom_line(aes(y = value, x = year, group = var),
    #              alpha = 0.1, size = 0.2) +
    #   facet_wrap(~common, scales = "free") +
    #   geom_line(data = Y.SR %>%
    #               mutate(year = 2002:2018) %>%
    #               gather(common, value, -year),
    #             aes(y = value,
    #                 x = year),
    #             color = "darkorange",
    #             size = 0.5) +
    #   labs(y = "Biomass density (kg km<sup>-2</sup>)") +
    #   theme_minimal() +
    #   theme(axis.title.x = element_blank(),
    #         axis.title.y = element_markdown())
    # 
    # ggsave(simulated_metapops,
    #        filename = here::here("figs/simulated_metapops.png"),
    #        bg = "white",
    #        dpi = 300,
    #        width = 6.5, height = 6)
    prog.bar=txtProgressBar(min=0, max=nrands, style=3)
    for (i in 1:length(unique(sim_out$var))){
      rand.mat <- 
        sim_out %>% 
        filter(var == unique(var)[i]) %>% 
        # dplyr::select(-var) %>% 
        group_by(grp) %>%
        mutate(year = rep(1:length(grp))) %>% 
        dplyr::select(-var) %>% 
        spread(.,grp,value) %>% 
        dplyr::select(-year)
      
      phi_SC_R_rands[i]=sync.fun(rand.mat)
      setTxtProgressBar(prog.bar, i)
    }
    
  } else {
    for (i in 1:nrands){
      rand.mat=cs_perm(Y.SR, ng = 1)
      phi_SC_R_rands[i]=sync.fun(rand.mat)
      setTxtProgressBar(prog.bar, i)
    }
    cat("\n")
  }

  # formatting results
  phi_SC_R_rands[nrands+1] <- obs$phi_SC_R
  # ses(nrands = length(phi_SC_R_rands))$plt
  phi_SC_R_pval_greater <- sum(phi_SC_R_rands >= obs$phi_SC_R)/(nrands+1)
  phi_SC_R_pval_lower <- sum(phi_SC_R_rands <= obs$phi_SC_R)/(nrands+1)
  phi_SC_R_lw95 <- quantile(phi_SC_R_rands, 0.025)
  phi_SC_R_up95 <- quantile(phi_SC_R_rands, 0.975)
  (phi_SC_R_SES <- (obs$phi_SC_R-mean(phi_SC_R_rands))/sd(phi_SC_R_rands))
  print(paste("phi_SC_R_SES:",phi_SC_R_SES))
  phi_SC_R_perc <- (obs$phi_SC_R - mean(phi_SC_R_rands))/obs$phi_SC_R*100
  
  # testing community-level spatial synchrony
  cat("Testing community-level spatial synchrony \n")
  phi_C_LR_rands <- numeric(length=nrands+1)*NA
  prog.bar=txtProgressBar(min=0, max=nrands, style=3)
  if (fit_gams){
    mods <- cs_perm2(Y.CL, ng = ncol(Y.CL), series = nrow(Y.CL),
                     fit_trend)
    sim_out <- NULL
    for (j in 1:nrow(mods)){
      # extract time series model for each patch
       m <- 
        mods %>% 
        ungroup() %>% 
        slice(j) %>% 
        .$sim
      
      if (!all(is.na(m[[1]]))){
        p <- m[[1]]$p
        mu <- m[[1]]$mu
        
        sim <-  replicate(n = nrands,
                          rTweedie(mu, p) + 
                            exp(as.numeric(simulate(m[[1]])))
        ) %>% 
          as.data.frame() %>% 
          gather(var, value) %>% 
          mutate(grp = j)
      } else {
        sim <- tibble(var = rep(paste0("V",1:nrands), each = nrow(mods$data[[1]])),
                      value = 0,
                      grp = j)
      }
       assign("sim_out", rbind(sim_out, sim))
       setTxtProgressBar(prog.bar, j)
    }
    # ggplot(sim_out %>%
    #         filter(grp %in% 1:10) %>% 
    #         group_by(grp, var) %>%
    #         mutate(year = 2002:2018)) +
    #          geom_line(aes(y = value, x = year, group = var),
    #                     alpha = 0.1, size = 0.2) +
    #          facet_wrap(~grp, scales = "free") 
    #          geom_line(data = Y.SR %>%
    #                      mutate(year = 2002:2018) %>%
    #                      gather(common, value, -year),
    #                    aes(y = value,
    #                        x = year),
    #                    color = "darkorange",
    #                    size = 0.5) +
    #          labs(y = "Biomass density (kg km<sup>-2</sup>)") +
    #          theme_minimal() +
    #          theme(axis.title.x = element_blank(),
    #                axis.title.y = element_markdown())
    prog.bar=txtProgressBar(min=0, max=nrands, style=3)
    for (i in 1:length(unique(sim_out$var))){
      rand.mat <- 
        sim_out %>% 
        filter(var == unique(var)[i]) %>% 
        # dplyr::select(-var) %>% 
        group_by(grp) %>%
        mutate(year = rep(1:length(grp))) %>% 
        dplyr::select(-var) %>% 
        spread(.,grp,value) %>% 
        dplyr::select(-year)
      
      phi_C_LR_rands[i]=sync.fun(rand.mat)
      setTxtProgressBar(prog.bar, i)
    }
    
  } else {
    prog.bar=txtProgressBar(min=0, max=nrands, style=3)
    for (i in 1:nrands){
      rand.mat=cs_perm(Y.CL, ng = 1)
      phi_C_LR_rands[i]=sync.fun(rand.mat)
      setTxtProgressBar(prog.bar, i)
    }
  }

  cat("\n")
  # formatting results

  # ses(phis = phi_C_LR_rands,
  #     nrands = length(phi_C_LR_rands),
  #     obs2 = obs$phi_C_LR)$plt
  phi_C_LR_rands[nrands+1] <- obs$phi_C_LR
  phi_C_LR_pval_greater <- sum(phi_C_LR_rands >= obs$phi_C_LR)/(nrands+1)
  phi_C_LR_pval_lower <- sum(phi_C_LR_rands <= obs$phi_C_LR)/(nrands+1)
  phi_C_LR_lw95 <- quantile(phi_C_LR_rands, 0.025)
  phi_C_LR_up95 <- quantile(phi_C_LR_rands, 0.975)
  phi_C_LR_SES <- (obs$phi_C_LR-mean(phi_C_LR_rands))/sd(phi_C_LR_rands)
  print(paste("phi_C_LR_SES:",phi_C_LR_SES))
  phi_C_LR_perc <- (obs$phi_C_LR - mean(phi_C_LR_rands))/obs$phi_C_LR*100
  
  # Overall test for:
  # average species-level spatial synchrony (phi_S,L->R)
  phi_S_LR_rands <- numeric(length=nrands+1)*NA
  # average local-scale species synchrony (phi_S->C,L)
  phi_SC_L_rands <- numeric(length=nrands+1)*NA
  # species-level spatial synchrony (phi_i,L->R)
  phi_i_LR_rands <- array(NA, dim=c(dim(Y)[2],2,nrands+1))
  # patch-level species synchrony (phi_S->C,k)
  phi_SC_k_rands <- array(NA, dim=c(s,2,nrands+1))
  cat("Overall test for average species-level spatial synchrony and average local-scale species synchrony \n")
  prog.bar <- txtProgressBar(min=0, max=nrands, style=3)
  if (fit_gams){
    SiteL.mods <- lapply(SiteL, function(x){ cs_perm2(x, ng = ncol(x), 
                                                      series = nrow(x),
                                                      fit_trend) })
    # SiteL.mods <- list(SiteL.mods[[1]],
    #                    SiteL.mods[[2]])
    sim_out <- NULL
    for (k in 1:length(SiteL.mods)){
      # for each species within each site...
      mods <- SiteL.mods[[k]] %>% ungroup()
      for (j in 1:nrow(mods)){
        # extract time series model for each species
        m <- 
          mods %>% 
          ungroup() %>% 
          slice(j) %>% 
          .$sim
        
        species <- mods %>% 
          ungroup() %>% 
          slice(j) %>% 
          pull(var)
        
        if (!all(is.na(m[[1]]))){
          p <- m[[1]]$p
          mu <- m[[1]]$mu
          
          sim <-  replicate(n = nrands,
                            rTweedie(mu, p) + 
                              exp(as.numeric(simulate(m[[1]])))
          ) %>% 
            as.data.frame() %>% 
            gather(sim, value) %>% 
            mutate(patch = k,
                   species,
                   grp = j)
        } else {
          sim <- tibble(sim = rep(paste0("V",1:nrands), each = nrow(mods$data[[1]])),
                        value = 0,
                        patch = k,
                        species,
                        grp = j)
        }

        
        assign("sim_out", rbind(sim_out, sim))
      }
      setTxtProgressBar(prog.bar, k)
    }
    
    # each facet is a different species in one of two patches. 
    # lines are different simulations (N = 1000)
    # ggplot(sim_out %>% 
    #          filter(sim == unique(sim)[i]) %>% 
    #          group_by(grp, sim, species, patch) %>%
    #          mutate(year = 2002:2018)) +
    #   geom_line(aes(y = value, x = year, group = sim)) +
    #   facet_wrap(patch~species)
    prog.bar <- txtProgressBar(min=0, max=nrands, style=3)
    for (i in 1:length(unique(sim_out$sim))){
      # each simulation gives all species time series at all locations
      rand.mat <-  sim_out %>% 
        filter(sim == unique(sim)[i]) %>%  
        group_by(patch, species) %>%
        mutate(year = rep(1:length(patch))) %>%
        ungroup() %>% 
        dplyr::select(-grp,-sim) %>% 
        # {. ->> vis_out} %>% 
        spread(.,species,value) %>% 
        dplyr::select(-1, -2)
    
      rand.res <- part_stab_comp(rand.mat, s = k, t)
      phi_S_LR_rands[i]=rand.res$part$phi_S_LR
      phi_SC_L_rands[i]=rand.res$part$phi_SC_L
      phi_i_LR_rands[,,i] <- as.matrix(rand.res$spe)
      phi_SC_k_rands[,,i] <- as.matrix(rand.res$patch)
      setTxtProgressBar(prog.bar, i)
    }
    
  } else {
    for (i in 1:nrands){
      # permute within-community time series
      SiteL.rand <- lapply(SiteL, function(x){ cs_perm(x, ng = 1) })
      # bind permutations into one big matrx (nrow == n patches * n years)
      Y.rand <- do.call("rbind", SiteL.rand)
      rand.res <- part_stab_comp(Y.rand, s, t)
      phi_S_LR_rands[i]=rand.res$part$phi_S_LR
      phi_SC_L_rands[i]=rand.res$part$phi_SC_L
      phi_i_LR_rands[,,i] <- as.matrix(rand.res$spe)
      phi_SC_k_rands[,,i] <- as.matrix(rand.res$patch)
      setTxtProgressBar(prog.bar, i)
    }
  }
  
  
  
  cat("\n")
  # formatting results
  phi_S_LR_rands[nrands+1] <- obs$phi_S_LR
  # ses(phis = phi_S_LR_rands,
  #     nrands = length(phi_S_LR_rands),
  #     obs2 = obs$phi_S_LR)$plt
  phi_S_LR_pval_greater <- sum(phi_S_LR_rands >= obs$phi_S_LR)/(nrands+1)
  phi_S_LR_pval_lower <- sum(phi_S_LR_rands <= obs$phi_S_LR)/(nrands+1)
  phi_S_LR_lw95 <- quantile(phi_S_LR_rands, 0.025)
  phi_S_LR_up95 <- quantile(phi_S_LR_rands, 0.975)
  phi_S_LR_SES <- (obs$phi_S_LR-mean(phi_S_LR_rands))/sd(phi_S_LR_rands)
  phi_S_LR_perc <- (obs$phi_S_LR-mean(phi_S_LR_rands))/obs$phi_S_LR*100
  
  # formatting results
  phi_SC_L_rands[nrands+1] <- obs$phi_SC_L
  phi_SC_L_pval_greater <- sum(phi_SC_L_rands >= obs$phi_SC_L)/(nrands+1)
  phi_SC_L_pval_lower <- sum(phi_SC_L_rands <= obs$phi_SC_L)/(nrands+1)
  phi_SC_L_lw95 <- quantile(phi_SC_L_rands, 0.025)
  phi_SC_L_up95 <- quantile(phi_SC_L_rands, 0.975)
  phi_SC_L_SES <- (obs$phi_SC_L-mean(phi_SC_L_rands))/sd(phi_SC_L_rands)
  phi_SC_L_perc <- (obs$phi_SC_L-mean(phi_SC_L_rands))/obs$phi_SC_L*100
  
  # formatting results for species-level spatial synchrony (phi_i,L->R) 
  phi_i_LR_rands[,,nrands+1] <- obs_spe
  res_spe <- cbind(obs_spe, data.frame(
    species=colnames(Y),
    pval_greater=apply(phi_i_LR_rands, 1:2, function(x){ sum(x >= x[nrands+1])/(nrands+1) })[,1],
    pval_lower=apply(phi_i_LR_rands, 1:2, function(x){ sum(x <= x[nrands+1])/(nrands+1) })[,1],
    lw95=apply(phi_i_LR_rands, 1:2, function(x){ quantile(x, 0.025) })[,1],
    up95=apply(phi_i_LR_rands, 1:2, function(x){ quantile(x, 0.975) })[,1],
    phi_i_LR_SES=(obs_spe[,1]-apply(phi_i_LR_rands, 1:2, mean)[,1])/apply(phi_i_LR_rands, 1:2, sd)[,1],
    phi_i_LR_perc=(obs_spe[,1]-apply(phi_i_LR_rands, 1:2, mean)[,1])/obs_spe[,1]*100))
  
  # formatting results for patch-level species synchrony (phi_S->C,k) 
  phi_SC_k_rands[,,nrands+1] <- obs_patch
  res_patch <- cbind(obs_patch, data.frame(
    patch=as.factor(1:s),
    pval_greater=apply(phi_SC_k_rands, 1:2, function(x){ sum(x >= x[nrands+1])/(nrands+1) })[,1],
    pval_lower=apply(phi_SC_k_rands, 1:2, function(x){ sum(x <= x[nrands+1])/(nrands+1) })[,1],
    lw95=apply(phi_SC_k_rands, 1:2, function(x){ quantile(x, 0.025) })[,1],
    up95=apply(phi_SC_k_rands, 1:2, function(x){ quantile(x, 0.975) })[,1],
    phi_SC_k_SES=(obs_patch[,1]-apply(phi_SC_k_rands, 1:2, mean)[,1])/apply(phi_SC_k_rands, 1:2, sd)[,1],
    phi_SC_k_perc=(obs_patch[,1]-apply(phi_SC_k_rands, 1:2, mean)[,1])/obs_patch[,1]*100))
  
  # summary result
  sync <- data.frame(sync=c("phi_SC_L", "phi_C_LR", "phi_S_LR", "phi_SC_R"),
                     pval_greater=c(phi_SC_L_pval_greater, phi_C_LR_pval_greater, phi_S_LR_pval_greater, phi_SC_R_pval_greater),
                     pval_lower=c(phi_SC_L_pval_lower, phi_C_LR_pval_lower, phi_S_LR_pval_lower, phi_SC_R_pval_lower),
                     lw95=c(phi_SC_L_lw95, phi_C_LR_lw95, phi_S_LR_lw95, phi_SC_R_lw95),
                     up95=c(phi_SC_L_up95, phi_C_LR_up95, phi_S_LR_up95, phi_SC_R_up95),
                     SES=c(phi_SC_L_SES, phi_C_LR_SES, phi_S_LR_SES, phi_SC_R_SES),
                     obs=c(obs$phi_SC_L, obs$phi_C_LR, obs$phi_S_LR, obs$phi_SC_R),
                     null = c(mean(phi_SC_L_rands), mean(phi_C_LR_rands), mean(phi_S_LR_rands), mean(phi_SC_R_rands)),
                     perc=c(phi_SC_L_perc, phi_C_LR_perc, phi_S_LR_perc, phi_SC_R_perc))
  # list of results
  res <- list(sync=sync, spe=res_spe, patch=res_patch)
  
  res
}
