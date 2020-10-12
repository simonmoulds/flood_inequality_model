gini = function(x) {
    if (length(which(!is.na(x))) < 2) {
        return(NA)
    }     
    if (min(x) < 0) {
        x = x - min(x)
    }
    x = x + 0.0000001
    x = sort(x)
    rank = order(x) # ascending order
    n = length(x)
    return(sum((2 * rank - n - 1) * x) / (n * sum(x)))
}

damage_function = function(W, xi_H, H, alpha_H, D) {
    F = 1 - exp(-(W + xi_H * H) / (alpha_H * D))
    F
}

compute_gini = function(x) {
    ns = dim(x)[1]
    nt = dim(x)[2]
    nc = dim(x)[3]    
    gini_index = apply(x, MARGIN=c(1,2), FUN=gini)
    dimnames(gini_index) = list(
        "scenario" = sprintf("s%d", 1:ns),
        "time" = 1:nt
    )
    gini_index
}

compute_combined_wealth = function(x) {
    ns = dim(x)[1]
    nt = dim(x)[2]
    nc = dim(x)[3]    
    wealth = apply(x, MARGIN=c(1,2), FUN=sum, na.rm=TRUE)
    wealth = t(wealth)
    dimnames(wealth) = list(
        "time" = 1:nt,
        "scenario" = sprintf("s%d", 1:ns)
    )
    wealth
}

get_water_level_ts = function(maxtime=50, dt=0.01, H_mult=1, theta=0.28) {
    t = seq(0, maxtime, dt)
    n_time = length(t)
    ## Inter-arrival time of flood events, plus flood magnitude
    t_mean = 1/dt
    w = rep(0, n_time)
    index = 0
    t0 = 0
    t1 = 0
    while (index < n_time) {
        x = runif(1)
        t1 = t0 + (t_mean * log(1/(1-x)))
        index = which.min(abs(t/dt - t1))
        y = runif(1)
        w[index] = H_mult * ((theta + 1) / theta) * (1 - (1 - y) ** theta)
        t0 = t1
    }
    w
}

get_default_parameter_set = function(ns, nc) {
    params = list(
        xi_H      = array(rep(c(0.5,0.5), each=ns), c(ns, nc)),
        gamma_E   = array(rep(c(0.005,0.005), each=ns), c(ns, nc)),
        epsilon_T = array(rep(c(1.1,1.1), each=ns), c(ns, nc)),
        keta_T    = array(rep(c(0.1,0.1), each=ns), c(ns, nc)),
        alpha_H   = array(rep(c(10,10), each=ns), c(ns, nc)),
        rho_E     = array(rep(c(1,1), each=ns), c(ns, nc)),
        lambda_P  = array(rep(c(2,2), each=ns), c(ns, nc)),
        phi_P   = array(rep(c(0.1,0.1), each=ns), c(ns, nc)),
        alpha_S   = array(rep(c(0.5,0.5), each=ns), c(ns, nc)),
        mu_S      = array(rep(c(0.25,0.25), each=ns), c(ns, nc)),
        tau       = array(0, c(ns, nc)),
        mu_S_prime= array(rep(c(0.5,0.5), each=ns), c(ns, nc))
    )
    params
}

get_default_initial_values = function(nc) {
    return(
        list(
            G=rep(0.01, nc),
            D=rep(1, nc),
            H=rep(0, nc),
            M=rep(0, nc)
        )
    )
}

plot_wealth_gini = function(wealth, gini_index, param_combinations, maxtime, dt) {
    
    Gini_cube = 
        as.tbl_cube(gini_index) %>% 
        as_tibble %>%
        mutate_if(sapply(., is.character), as.factor) %>%
        left_join(param_combinations, by="scenario")    
    
    Wealth_cube =
        as.tbl_cube(wealth) %>% 
        as_tibble %>%
        mutate_if(sapply(., is.character), as.factor) %>%
        left_join(param_combinations, by="scenario")

    data_cube = Gini_cube %>% left_join(Wealth_cube)
    
    lambda_P_vals =
        unique(param_combinations$lambda_P) %>%
        sort %>%
        format(digits=4, drop0trailing = TRUE)
    alpha_H_vals =
        unique(param_combinations$alpha_H) %>%
        sort %>%
        format(digits=4, drop0trailing=TRUE)
    data_cube$lambda_P = factor(
        data_cube$lambda_P,
        labels = sapply(lambda_P_vals, FUN=function(x) bquote(lambda[P]==.(x)))
    )
    data_cube$alpha_H = factor(
        data_cube$alpha_H,
        labels = sapply(alpha_H_vals, FUN=function(x) bquote(alpha[H]==.(x)))
    )
    data_cube$community = factor(data_cube$community)
    data_cube$wealth = log10(data_cube$wealth)
    dummy_wealth = data_cube %>%
        dplyr::select(-gini_index) %>%
        spread(key=community, value=wealth) %>%
        mutate(`3`=0) %>%
        gather(community, wealth, `1`, `2`, `3`)
    
    Gini_cube_limits =
        data_cube %>%
        dplyr::select(-wealth) %>% 
        spread(key=runs, value=gini_index) %>%
        mutate(
            gini_max=apply(.[grep("^r", names(.))], 1, max, na.rm=TRUE),
            gini_min=apply(.[grep("^r", names(.))], 1, min, na.rm=TRUE)
        ) %>%
        mutate(
            gini_max = gini_max * 20 - 5, # scale
            gini_min = gini_min * 20 - 5  # scale
        ) %>%
        mutate(
            gini_max = ifelse(is.finite(gini_max), gini_max, NA),
            gini_min = ifelse(is.finite(gini_min), gini_min, NA)
        )
    
    Gini_cube_limits$community %<>% as.character

    ## NEW:
    Gini_cube_limits %<>% dplyr::select(!starts_with('r'))
    
    data_cube %<>%
        left_join(Gini_cube_limits)

    data_cube %<>% mutate(gini_index = gini_index * 20 - 5)
    ## NEW:
    ## data_cube %<>%
    ##     right_join(dummy_wealth)
    ## data_cube %<>% select(-gini_max, -gini_min)

    p =
        ggplot(
            data=data_cube,
            aes(
                x=time,
                y=wealth,
                ## ymin=gini_min,
                ## ymax=gini_max,
                group=interaction(runs,community),
                color=community
            )
        ) +
        ## geom_line(size=0.05, alpha=.75) + 
        geom_line(
            aes(
                x=time,
                y=gini_index,
                color="black"
            ),
            size=0.01, alpha=.5
        ) +
        geom_line(size=0.025, alpha=1) + #.75) + 
        ylab(expression(log[10]*"(G)")) +
        xlab("Time") + 
        scale_x_continuous(
            breaks=seq(0,maxtime/dt,maxtime/dt),
            labels=c("0","50")
        ) +
        scale_y_continuous(
            breaks=seq(-5, 5, 10),
            limits=c(-5,5),
            sec.axis=sec_axis(
                ~(.+5)/20,
                breaks=seq(0,0.5,0.5),
                name="Gini coefficient"
            )
        ) +
        scale_color_manual(
            name=element_blank(),
            labels=c("Planned","Unplanned","Gini coefficient"),
            values=c("magenta","cyan3","gray40")
        ) +
        ## scale_fill_discrete(guide="none") + 
        facet_grid(lambda_P ~ alpha_H, labeller=label_parsed) +
        theme(
            aspect.ratio=1,
            axis.line = element_line(color="black", size=0.25),
            axis.text.x=element_text(size=rel(0.8)),
            axis.text.y=element_text(size=rel(0.8)),
            strip.text.x=element_text(size=rel(0.8)),
            strip.text.y=element_text(size=rel(0.8)),
            axis.title.x=element_text(size=rel(0.8)),
            axis.title.y=element_text(size=rel(0.8)),
            panel.background = element_rect(
                fill="transparent",
                colour="black",
                size=0.25,
                linetype=1),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            strip.background=element_blank(),
            legend.position=c(0.5,-0.15)
        ) +        
        guides(
            color=guide_legend(
                override.aes=list(
                    color=c("magenta","cyan3","gray40"),
                    size=0.5,
                    fill=c("white","white","white"),
                    alpha=c(1,1,1)), ncol=3, byrow=TRUE
            )
        )
    p
}

plot_wealth_gini2 = function(wealth, gini_index, param_combinations, maxtime, dt) {
    ## Same as plot_wealth_gini, but adjusted so the mu_S_prime varies
    ## in vertical direction

    Gini_cube = 
        as.tbl_cube(gini_index) %>% 
        as_tibble %>%
        mutate_if(sapply(., is.character), as.factor) %>%
        left_join(param_combinations, by="scenario")    
    
    Wealth_cube =
        as.tbl_cube(wealth) %>% 
        as_tibble %>%
        mutate_if(sapply(., is.character), as.factor) %>%
        left_join(param_combinations, by="scenario")

    data_cube = Gini_cube %>% left_join(Wealth_cube)
    
    mu_S_prime_vals =
        unique(param_combinations$mu_S_prime) %>%
        sort %>%
        format(digits=4, drop0trailing = TRUE)
    alpha_H_vals =
        unique(param_combinations$alpha_H) %>%
        sort %>%
        format(digits=4, drop0trailing=TRUE)
    data_cube$mu_S_prime = factor(
        data_cube$mu_S_prime,
        labels = sapply(mu_S_prime_vals, FUN=function(x) bquote(paste(plain(mu[S]), plain("'")==.(x))))
        ## labels = sapply(mu_S_prime_vals, FUN=function(x) bquote(lambda[P]==.(x)))
    )
    data_cube$alpha_H = factor(
        data_cube$alpha_H,
        labels = sapply(alpha_H_vals, FUN=function(x) bquote(alpha[H]==.(x)))
    )
    data_cube$community = factor(data_cube$community)
    data_cube$wealth = log10(data_cube$wealth)
    ## dummy_wealth = data_cube %>% dplyr::select(-gini_index) %>% spread(key=community, value=wealth) %>% mutate(`3`=0) %>% gather(community, wealth, `1`, `2`, `3`)
    
    Gini_cube_limits =
        data_cube %>%
        dplyr::select(-wealth) %>% 
        spread(key=runs, value=gini_index) %>%
        mutate(
            gini_max=apply(.[grep("^r", names(.))], 1, max, na.rm=TRUE),
            gini_min=apply(.[grep("^r", names(.))], 1, min, na.rm=TRUE)
        ) %>%
        mutate(
            gini_max = gini_max * 20 - 5, # scale
            gini_min = gini_min * 20 - 5  # scale
        ) %>%
        mutate(
            gini_max = ifelse(is.finite(gini_max), gini_max, NA),
            gini_min = ifelse(is.finite(gini_min), gini_min, NA)
        )
    
    Gini_cube_limits$community %<>% as.character    
    Gini_cube_limits %<>% dplyr::select(!starts_with('r'))   
    data_cube %<>%
        left_join(Gini_cube_limits)
    data_cube %<>% mutate(gini_index = gini_index * 20 - 5)

    p =
        ggplot(
            data=data_cube,
            aes(
                x=time,
                y=wealth,
                ## ymin=gini_min,
                ## ymax=gini_max,
                group=interaction(runs,community),
                color=community
            )
        ) +
        ## geom_line(size=0.05, alpha=.75) + 
        geom_line(
            aes(
                x=time,
                y=gini_index,
                color="black"
            ),
            size=0.01, alpha=.5
        ) +
        geom_line(size=0.025, alpha=1) + #.75) + 
        ylab(expression(log[10]*"(G)")) +
        xlab("Time") + 
        scale_x_continuous(
            breaks=seq(0,maxtime/dt,maxtime/dt),
            labels=c("0","50")
        ) +
        scale_y_continuous(
            breaks=seq(-5, 5, 10),
            limits=c(-5,5),
            sec.axis=sec_axis(
                ~(.+5)/20,
                breaks=seq(0,0.5,0.5),
                name="Gini coefficient"
            )
        ) +
        scale_color_manual(
            name=element_blank(),
            labels=c("Planned","Unplanned","Gini coefficient"),
            values=c("magenta","cyan3","gray40")
        ) +
        ## scale_fill_discrete(guide="none") + 
        facet_grid(mu_S_prime ~ alpha_H, labeller=label_parsed) +
        ## facet_grid(lambda_P ~ alpha_H, labeller=label_parsed) +
        theme(
            aspect.ratio=1,
            axis.line = element_line(color="black", size=0.25),
            axis.text.x=element_text(size=rel(0.8)),
            axis.text.y=element_text(size=rel(0.8)),
            strip.text.x=element_text(size=rel(0.8)),
            strip.text.y=element_text(size=rel(0.8)),
            axis.title.x=element_text(size=rel(0.8)),
            axis.title.y=element_text(size=rel(0.8)),
            panel.background = element_rect(
                fill="transparent",
                colour="black",
                size=0.25,
                linetype=1),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            strip.background=element_blank(),
            legend.position=c(0.5,-0.15)
        ) +        
        guides(
            color=guide_legend(
                override.aes=list(
                    color=c("magenta","cyan3","gray40"),
                    size=0.5,
                    fill=c("white","white","white"),
                    alpha=c(1,1,1)), ncol=3, byrow=TRUE
            )
        )
    p
}

## constrain_protection_increase = function(R0, F, G, H, gamma_E, keta_T, check_incentive=TRUE, check_affordable=TRUE, allow_suboptimal=FALSE) {
##     ## Function to balance the desired flood protection increase
##     ## against the community socioeconomics.
##     ##
##     ## Args:
##     ## =====
##     ## R0 : Desired increase in flood protection
##     ## F  : Flood damage experienced by the community
##     ## G  : Community wealth
##     ## gamma_E : Unit cost of flood protection increase
##     ## check_incentive :
##     ## check_affordable :
##     ## allow_suboptimal : if TRUE, allow community to build flood
##     ##   protection which would not protect against flood damage F
##     ##
##     ## Return:
##     ## =======
##     ## Increase in flood protection considering community
##     ##   socioeconomics
##     ##    
##     R = R0

##     if (is.na(G)) {
##         return(0)
##     }    
    
##     ## 1. Is there an economic incentive to install the levee?
##     if (check_incentive) {
##         if (!((F * G) >= (gamma_E * R * sqrt(G)))) {
##             R = 0                
##         }
##     }
    
##     ## 2. Is there enough wealth to install the levee?
##     if (check_affordable) {
##         if (!((G - F * G) >= (gamma_E * R * sqrt(G)))) {
##         ## if (!(T >= (gamma_E * R * sqrt(G)))) {
##             R = 0            
##             if (allow_suboptimal) {
##                 ## determine the maximum levee height that can be installed
##                 R_max = (G - F * G) / (gamma_E * sqrt(G))
##                 R = R_max
##             }        
##         }
##     }

##     ## ## 3. Does the height of the levee exceed an allowable amount?
##     ## R_max = H_max - H + keta_T * H
##     ## if (R > R_max) {
##     ##     R = R_max
##     ## }
##     R
## }

compute_protection_increase = function(F, G, H, W, epsilon_T, xi_H, gamma_E, keta_T, check_incentive, check_affordable, allow_suboptimal=FALSE) {
    ## Function to compute flood protection increase.
    ##
    ## Args:
    ## =====
    ## F : Flood damage
    ## W : Water level
    ## H : Current protection level    
    ## xi_H : TODO
    ##
    ## Return:
    ## =======
    ## Increase in flood protection
    ##
    if (F > 0) {
        R = epsilon_T * (W + xi_H * H - H)
    } else {
        R = 0
    }
    if (is.na(G)) {
        return(0)
    }        
    ## 1. Is there an economic incentive to install the levee?
    if (check_incentive) {
        if (!((F * G) >= (gamma_E * R * sqrt(max(G, 0))))) {
            R = 0                
        }
    }    
    ## 2. Is there enough wealth to install the levee?
    if (check_affordable) {
        if (!((G - F * G) >= (gamma_E * R * sqrt(max(G, 0))))) {
        ## if (!(T >= (gamma_E * R * sqrt(G)))) {
            R = 0            
            if (allow_suboptimal) {
                ## determine the maximum levee height that can be installed
                R_max = (G - F * G) / (gamma_E * sqrt(G))
                R = R_max
            }        
        }
    }
    ## ## 3. Does the height of the levee exceed an allowable amount?
    ## R_max = H_max - H + keta_T * H
    ## if (R > R_max) {
    ##     R = R_max
    ## }
    R    
}

model_ineq = function(w, init, params, n_community, n_scenario) {    

    n_time = length(w)
    
    ## set dimensions
    dims = list(
        "scenario" = sprintf("s%d", 1:n_scenario),
        "time" = 1:n_time,
        "community" = sprintf("c%d", 1:n_community)
    )
    ndims = sapply(dims, length)
    
    ## prepare output variables
    F = array(NA, ndims) %>% `dimnames<-`(dims)
    R = array(NA, ndims) %>% `dimnames<-`(dims)
    R0 = array(NA, ndims) %>% `dimnames<-`(dims)
    S = array(NA, ndims) %>% `dimnames<-`(dims)
    G = array(NA, ndims) %>% `dimnames<-`(dims)
    D = array(NA, ndims) %>% `dimnames<-`(dims)
    H = array(NA, ndims) %>% `dimnames<-`(dims)
    M = array(NA, ndims) %>% `dimnames<-`(dims)
    G[,1,] = rep(init$G, each=n_scenario)
    D[,1,] = rep(init$D, each=n_scenario)
    H[,1,] = rep(init$H, each=n_scenario)
    M[,1,] = rep(init$M, each=n_scenario)

    ## TEST
    T = array(NA, ndims) %>% `dimnames<-`(dims)
    T[,1,] = 0
    Mp = array(NA, ndims) %>% `dimnames<-`(dims)
    Mp[,1,] = rep(init$M, each=n_scenario)
    
    ## run simulation
    for (i in 1:n_scenario) {

        failure = rep(FALSE, n_community)            
        
        for (j in 1:(n_time - 1)) {

            ## value of Dirac comb - does flooding occur?
            if (w[j] > 0) {
                Delta = 1 / dt
            } else {
                Delta = 0
            }

            ## First compute damage and initial values of
            ## flood protection increase and shock
            for (k in 1:n_community) {

                if (!failure[k]) {
                    
                    ## damage function: proportion of settlement
                    ## that is damaged
                    if ((w[j] + params$xi_H[i,k] * H[i,j,k]) > H[i,j,k]) {
                        F[i,j,k] = 1 - exp(-(w[j] + params$xi_H[i,k] * H[i,j,k]) / (params$alpha_H[i,k] * D[i,j,k]))
                    } else {
                        F[i,j,k] = 0
                    }

                    ## raising of levees
                    R[i,j,k] = compute_protection_increase(
                        F[i,j,k],
                        G[i,j,k],
                        H[i,j,k],
                        w[j],
                        params$epsilon_T[i,k],
                        params$xi_H[i,k],
                        params$gamma_E[i,k],
                        params$keta_T[i,k],
                        check_incentive=TRUE,
                        check_affordable=TRUE
                    )

                    ## shock magnitude
                    if (R[i,j,k] > 0) {
                        S[i,j,k] = F[i,j,k] * params$alpha_S[i,k]
                    } else {
                        S[i,j,k] = F[i,j,k]
                    }                    
                }                                
            }
            
            ## now compute shock based on flooding experienced by
            ## other communities, and redistribution of wealth
            for (k in 1:n_community) {
                if (!failure[k]) {                    

                    ## Get shock and protection increase
                    Sk = S[i,j,]
                    Rk = R[i,j,]
                    ## set the shock level of the current community k
                    ## to zero, so that we're only considering the
                    ## shock experienced by other communities
                    Sk[k] = 0
                    Sk[failure] = 0 # don't transfer wealth to failed communities
                    ## calculate the change in awareness of flood events
                    ## outside of the current community
                    ## - is max(Sk) the best way to summarise?
                    ## - in an ABM, this could be scaled by a nb function
                    Sp = max(Sk) # is max(Sk) the best way to summarise?
                    ## update memory of flood events which occur
                    ## in other communities
                    ## ***mu_S should be replaced***
                    dMp_dt = Delta * Sp - params$mu_S_prime[i,k] * Mp[i,j,k]
                    ## Mp[i,j+1,k] = Mp[i,j,k] + dMp_dt * dt
                    Mp[i,j+1,k] = max(0, Mp[i,j,k] + dMp_dt * dt)
                    
                    ## Compute redistribution, then check its affordable
                    ## by assuming that a community will first install
                    ## flood protection itself before helping others
                    ## ***need to check that this is robust when dealing with failed communities***
                    Gk = G[i,j,]
                    T[i,j,k] = (mean(Gk, na.rm=TRUE) - Gk[k]) * params$tau[i,k] * Mp[i,j+1,k]                    
                    ## initial change in wealth, without redistribution
                    ## N.B. we're assuming here that a community will
                    ## first install flood protection itself before
                    ## helping others; relax in some scenarios?
                    dG_dt_init = (
                        params$rho_E[i,k] * (1 - D[i,j,k]) * G[i,j,k]
                        - Delta * (
                            F[i,j,k]
                            * G[i,j,k]
                            + params$gamma_E[i,k]
                            * R[i,j,k]
                            * sqrt(G[i,j,k])
                        )
                    )
                    G_avail = G[i,j,k] + dG_dt_init * dt
                    T[i,j,k] = min(G_avail, T[i,j,k])
                } else {
                    T[i,j,k] = NA                    
                }                
            }

            ## scale T to account for the fact that the
            ## target redistribution communities with a net loss
            ## may have been scaled to account for affordability
            Tk = T[i,j,]            
            sf = abs(sum(Tk[Tk <= 0], na.rm=TRUE) / sum(Tk[Tk > 0], na.rm=TRUE))
            Tk[!is.na(Tk) & Tk > 0] %<>% `*`(sf)
            T[i,j,] = Tk
            
            for (k in 1:n_community) {
                if (!failure[k]) {
                    
                    ## after calculating the redistribution of wealth,
                    ## recalculate protection increase
                    R[i,j,k] = compute_protection_increase(
                        F[i,j,k],
                        G[i,j,k] + T[i,j,k],
                        H[i,j,k],
                        w[j],
                        params$epsilon_T[i,k],
                        params$xi_H[i,k],
                        params$gamma_E[i,k],
                        params$keta_T[i,k],
                        check_incentive=TRUE, # TODO: check this?
                        check_affordable=TRUE
                    )                    
                    ## shock magnitude                    
                    if (R[i,j,k] > 0) {
                        S[i,j,k] = F[i,j,k] * params$alpha_S[i,k]
                    } else {
                        S[i,j,k] = F[i,j,k]
                    }

                    ## ###################### ##
                    ## update state variables ##
                    ## ###################### ##

                    ## 1 - change in wealth
                    dG_dt = (
                        params$rho_E[i,k] * (1 - D[i,j,k]) * G[i,j,k]
                        + T[i,j,k]
                        - Delta * (
                            F[i,j,k]
                            * G[i,j,k]
                            + params$gamma_E[i,k]
                            * R[i,j,k]
                            * sqrt(G[i,j,k])
                        )
                    )                    
                    G[i,j+1,k] = max(0, G[i,j,k] + dG_dt * dt)
                    
                    ## 2 - change in distance from river
                    dD_dt = (
                        (M[i,j,k] - D[i,j,k] / params$lambda_P[i,k])
                        * (params$phi_P[i,k] / sqrt(G[i,j,k])
                        )
                    )                    
                    D[i,j+1,k] = max(0, D[i,j,k] + dD_dt * dt)

                    ## change in levee height
                    dH_dt = (
                        Delta * R[i,j,k] - params$keta_T[i,k] * H[i,j,k]
                    )                    
                    H[i,j+1,k] = max(0, H[i,j,k] + dH_dt * dt)

                    ## change in community memory
                    dM_dt = (
                        Delta * S[i,j,k] - params$mu_S[i,k] * M[i,j,k]
                    )                    
                    M[i,j+1,k] = max(0, M[i,j,k] + dM_dt * dt)

                    ## G[i,j+1,k] = max(0, G[i,j,k] + dG_dt * dt)
                    ## D[i,j+1,k] = max(0, D[i,j,k] + dD_dt * dt)
                    ## H[i,j+1,k] = max(0, H[i,j,k] + dH_dt * dt)
                    ## M[i,j+1,k] = max(0, M[i,j,k] + dM_dt * dt)

                    ## check failure
                    if (G[i,j+1,k] <= 0.00001) {
                        failure[k] = TRUE
                    }
                }
            }
        }
    }
    return(list(F=F, R=R, S=S, G=G, D=D, H=H, M=M, T=T))
}

## make_plot = function(G, H, D, F, labels) {
##     G_cube = as.tbl_cube(G) %>% as_tibble %>% mutate_if(sapply(., is.character), as.factor)
##     levels(G_cube$community) = c("Planned","Unplanned")
##     H_cube = as.tbl_cube(H) %>% as_tibble %>% mutate_if(sapply(., is.character), as.factor)
##     levels(H_cube$community) = c("Planned","Unplanned")
##     D_cube = as.tbl_cube(D) %>% as_tibble %>% mutate_if(sapply(., is.character), as.factor)
##     levels(D_cube$community) = c("Planned","Unplanned")
##     F_cube = as.tbl_cube(F) %>% as_tibble %>% mutate_if(sapply(., is.character), as.factor)
##     levels(F_cube$community) = c("Planned","Unplanned")
##     F_cube[F_cube == 0] = NA
##     p1 =
##         ggplot(data=G_cube, aes(x=time, y=log10(G), color=scenario)) +
##         facet_wrap(.~community) + 
##         geom_line() +
##         ylab("Wealth") +
##         xlab(element_blank()) +
##         ylim(c(-5,5)) +
##         scale_x_continuous(breaks=seq(0,5000,1000), labels=c("0","10","20","30","40","50")) +
##         theme(
##             strip.background=element_blank(),
##             plot.margin = margin(t=0,b=0),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             legend.position="none"
##             ## legend.position="bottom"
##         )
##     p2 =
##         ggplot(data=H_cube, aes(x=time, y=H, color=scenario)) +
##         facet_wrap(.~community) + 
##         ylab("Technology") +
##         xlab(element_blank()) +
##         ylim(c(0,5)) + 
##         scale_x_continuous(breaks=seq(0,5000,1000), labels=c("0","10","20","30","40","50")) +
##         geom_line() +
##         theme(
##             ## plot.margin = margin(t=0,b=0),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             strip.background=element_blank(),
##             strip.text.x=element_blank(),
##             legend.position="none"
##             ## strip.text.x=element_text(color="transparent")
##         )           
##     p3 =
##         ggplot(data=D_cube, aes(x=time, y=D, color=scenario)) +
##         facet_wrap(.~community) + 
##         ylab("Distance") +
##         xlab(element_blank()) +
##         ylim(c(0,2.5)) + 
##         scale_x_continuous(breaks=seq(0,5000,1000), labels=c("0","10","20","30","40","50")) +
##         geom_line() +
##         theme(
##             ## plot.margin = margin(t=0,b=0),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             strip.background=element_blank(),
##             strip.text.x=element_blank(),
##             legend.position="none"
##             ## strip.text.x=element_text(color="transparent")
##         )           
    
##     p4 =
##         ggplot(data=F_cube, aes(x=time, y=F, color=scenario)) +
##         facet_wrap(.~community) + 
##         ylab("Damage") +
##         xlab("Time") +
##         ## xlab(element_blank()) +
##         ylim(c(0,1)) + 
##         scale_x_continuous(breaks=seq(0,5000,1000), labels=c("0","10","20","30","40","50")) +
##         geom_point() +
##         geom_line(size=0.4, alpha=0) + 
##         scale_color_discrete(
##             name=element_blank(),
##             labels=labels,
##             ## labels=c(expression(paste(tau, "=0.")),
##             ##          expression(paste(tau, "=0.05")),
##             ##          expression(paste(tau, "=0.5")),
##             ##          expression(paste(tau, "=1."))
##             ##          )
##         ) +
##         theme(
##             ## plot.margin = margin(t=0,b=0),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             strip.background=element_blank(),
##             strip.text.x=element_blank(),
##             legend.position="bottom"
##             ## strip.text.x=element_text(color="transparent")
##         ) + guides(color = guide_legend(override.aes = list(alpha = 1), ncol=2, byrow=TRUE))
    
##     ## legend = get_legend(p1)
##     ## pcol = plot_grid(
##     ##     p1 + theme(legend.position="none"),
##     ##     p2 + theme(legend.position="none"),
##     ##     p3 + theme(legend.position="none"),
##     ##     p4 + theme(legend.position="none"),
##     ##     align="v",
##     ##     axis="bt",
##     ##     ## rel_heights=c(1.5,1,1,1),
##     ##     nrow=4
##     ## )
##     ## ## p = plot_grid(pcol,legend,ncol=1,rel_heights=c(1,.1)) +
##     ## ##     draw_label("Hello, world", x=0.33, y=0.9)
##     ## p1=p1 + theme(legend.position="none")
##     ## p2=p2 + theme(legend.position="none")
##     ## p3=p3 + theme(legend.position="none")
##     ## p4=p4 ##+ theme(legend.position="none")
##     ## p = ggarrange(p1,p2,p3,p4,nrow=4)
##     ## p = list(p1,p2,p3,p4)
##     pl <- lapply(list(p1,p2,p3,p4), ggplotGrob)
##     pl <- Reduce(gridExtra::gtable_rbind, pl)
##     ## grid::grid.draw(pl)
##     return(pl)
## }

## model_ineq = function(w, init, params, n_community, n_scenario) {    

##     n_time = length(w)
    
##     ## set dimensions
##     dims = list(
##         "scenario" = sprintf("s%d", 1:n_scenario),
##         "time" = 1:n_time,
##         "community" = sprintf("c%d", 1:n_community)
##     )
##     ndims = sapply(dims, length)
    
##     ## prepare output variables
##     F = array(NA, ndims) %>% `dimnames<-`(dims)
##     R = array(NA, ndims) %>% `dimnames<-`(dims)
##     R0 = array(NA, ndims) %>% `dimnames<-`(dims)
##     S = array(NA, ndims) %>% `dimnames<-`(dims)
##     G = array(NA, ndims) %>% `dimnames<-`(dims)
##     D = array(NA, ndims) %>% `dimnames<-`(dims)
##     H = array(NA, ndims) %>% `dimnames<-`(dims)
##     M = array(NA, ndims) %>% `dimnames<-`(dims)
##     G[,1,] = rep(init$G, each=n_scenario)
##     D[,1,] = rep(init$D, each=n_scenario)
##     H[,1,] = rep(init$H, each=n_scenario)
##     M[,1,] = rep(init$M, each=n_scenario)
    
##     ## run simulation
##     for (i in 1:n_scenario) {

##         failure = rep(FALSE, n_community)            
        
##         for (j in 1:(n_time - 1)) {
##             for (k in 1:n_community) {

##                 if (!failure[k]) {
                    
##                     ## damage function: proportion of settlement
##                     ## that is damaged
##                     if ((w[j] + params$xi_H[i,k] * H[i,j,k]) > H[i,j,k]) {
##                         F[i,j,k] = 1 - exp(-(w[j] + params$xi_H[i,k] * H[i,j,k]) / (params$alpha_H[i,k] * D[i,j,k]))
##                     } else {
##                         F[i,j,k] = 0
##                     }

##                     ## raising of levees: first calculate the theoretical
##                     ## amount by which the community would want to raise
##                     ## the levee, then constrain based on wealth/incentive
##                     R0[i,j,k] = compute_protection_increase(
##                         F[i,j,k],
##                         w[j],
##                         H[i,j,k],
##                         params$epsilon_T[i,k],
##                         params$xi_H[i,k]
##                     )
##                     R[i,j,k] = constrain_protection_increase(
##                         R0[i,j,k],
##                         F[i,j,k],
##                         G[i,j,k],
##                         H[i,j,k],
##                         ## H_max=params$H_max[i,k],
##                         params$gamma_E[i,k],
##                         params$keta_T[i,k],
##                         check_incentive=TRUE,
##                         check_affordable=TRUE
##                     )

##                     ## shock magnitude
##                     if (R[i,j,k] > 0) {
##                         S[i,j,k] = F[i,j,k] * params$alpha_S[i,k]
##                     } else {
##                         S[i,j,k] = F[i,j,k]
##                     }

##                     ## value of Dirac comb
##                     if (w[j] > 0) {
##                         Delta = 1 / dt
##                     } else {
##                         Delta = 0
##                     }

##                     ## compute system dynamics without redistribution
##                     income = params$rho_E[i,] * (1 - D[i,j,]) * G[i,j,]
##                     avg_income = mean(income, na.rm=TRUE)

##                     ## Hoover Index:
##                     redistr = rep(0, n_community)
##                     if (avg_income > 0) {
##                         HI = 0.5 * sum(abs(income - avg_income), na.rm=TRUE) / sum(income, na.rm=TRUE)
##                         if (HI > 0) {
##                             ## redistribution:
##                             rel = (income - mean(income, na.rm=TRUE)) / sum(abs(income - mean(income, na.rm=TRUE)), na.rm=TRUE) * 2
##                             redistr = sum(income) * HI * params$tau[i,k] * (income - mean(income)) / sum(abs(income - mean(income))) * -2
##                         }                
##                     }            

##                     ## change in wealth
##                     dG_dt = params$rho_E[i,k] * (1 - D[i,j,k]) * G[i,j,k] + redistr[k] - Delta * (F[i,j,k] * G[i,j,k] + params$gamma_E[i,k] * R[i,j,k] * sqrt(G[i,j,k]))
##                     ## dG_dt = params$rho_E[i,k] * (1 - D[i,j,k] ** (1/params$s[i,k])) * G[i,j,k] + redistr[k] - Delta * (F[i,j,k] * G[i,j,k] + params$gamma_E[i,k] * R[i,j,k] * sqrt(G[i,j,k]))

##                     ## change in distance from river
##                     dD_dt = (M[i,j,k] - D[i,j,k] / params$lambda_P[i,k]) * (params$phi_P[i,k] / sqrt(G[i,j,k]))

##                     ## change in levee height
##                     dH_dt = Delta * R[i,j,k] - params$keta_T[i,k] * H[i,j,k]

##                     ## change in community memory
##                     dM_dt = Delta * S[i,j,k] - params$mu_S[i,k] * M[i,j,k]

##                     ## update state variables
##                     G[i,j+1,k] = max(0, G[i,j,k] + dG_dt * dt)
##                     D[i,j+1,k] = max(0, D[i,j,k] + dD_dt * dt)
##                     H[i,j+1,k] = max(0, H[i,j,k] + dH_dt * dt)
##                     M[i,j+1,k] = max(0, M[i,j,k] + dM_dt * dt)

##                     ## check failure
##                     if (G[i,j+1,k] <= 0.00001) {
##                         failure[k] = TRUE
##                     }
##                 }                
##             }
##         }
##     }
##     return(list(F=F, R=R, S=S, G=G, D=D, H=H, M=M))
## }

## add_margin = function(...)  grid::convertUnit(theme_get()[["plot.margin"]] + margin(...), "cm")

## plot_gini = function(gini_index, param_combinations) {    
    
##     Gini_cube =
##         as.tbl_cube(gini_index) %>% 
##         as_tibble %>%
##         mutate_if(sapply(., is.character), as.factor) %>%
##         left_join(param_combinations, by="scenario")    
##     lambda_P_vals =
##         unique(param_combinations$lambda_P) %>%
##         sort %>%
##         format(digits=4, drop0trailing = TRUE)
##     alpha_H_vals =
##         unique(param_combinations$alpha_H) %>%
##         sort %>%
##         format(digits=4, drop0trailing=TRUE)
##     Gini_cube$lambda_P = factor(
##         Gini_cube$lambda_P,
##         labels = sapply(lambda_P_vals, FUN=function(x) bquote(lambda[P]==.(x)))
##     )
##     Gini_cube$alpha_H = factor(
##         Gini_cube$alpha_H,
##         labels = sapply(alpha_H_vals, FUN=function(x) bquote(alpha[H]==.(x)))
##     )    

##     p =
##         ggplot(data=Gini_cube, aes(x=time, y=gini_index, group=runs)) +
##         geom_line(size=0.05) +
##         ylab("Gini coefficient") +
##         xlab("Time") +
##         scale_x_continuous(
##             breaks=seq(0,5000,5000),
##             labels=c("0","50")
##         ) +
##         scale_y_continuous(
##             breaks=seq(0,0.5,0.5),
##             limits=c(0,0.5)## ,
##             ## sec.axis=sec_axis(~./5, name="Test")
##         ) +
##         facet_grid(lambda_P ~ alpha_H, labeller=label_parsed) +
##         theme(
##             aspect.ratio=1,
##             axis.text.x=element_text(size=rel(0.8)),
##             axis.text.y=element_text(size=rel(0.8)),
##             strip.text.x=element_text(size=rel(0.8)),
##             strip.text.y=element_text(size=rel(0.8)),
##             axis.title.x=element_text(size=rel(0.8)),
##             axis.title.y=element_text(size=rel(0.8)),
##             panel.background = element_rect(
##                 fill="transparent",
##                 colour="black",
##                 size=0.25,
##                 linetype=1),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             strip.background=element_blank(),
##             legend.position="bottom"
##         ) ## +
##         ## theme(plot.margin=add_margin(t=0, unit="cm"))
##     p
## }

## plot_wealth = function(wealth, param_combinations) {    
    
##     Wealth_cube =
##         as.tbl_cube(wealth) %>% 
##         as_tibble %>%
##         mutate_if(sapply(., is.character), as.factor) %>%
##         left_join(param_combinations, by="scenario")    
##     lambda_P_vals =
##         unique(param_combinations$lambda_P) %>%
##         sort %>%
##         format(digits=4, drop0trailing = TRUE)
##     alpha_H_vals =
##         unique(param_combinations$alpha_H) %>%
##         sort %>%
##         format(digits=4, drop0trailing=TRUE)
##     Wealth_cube$lambda_P = factor(
##         Wealth_cube$lambda_P,
##         labels = sapply(lambda_P_vals, FUN=function(x) bquote(lambda[P]==.(x)))
##     )
##     Wealth_cube$alpha_H = factor(
##         Wealth_cube$alpha_H,
##         labels = sapply(alpha_H_vals, FUN=function(x) bquote(alpha[H]==.(x)))
##     )
##     Wealth_cube$community = factor(Wealth_cube$community)
##     Wealth_cube$wealth = log10(Wealth_cube$wealth)
##     p =
##         ggplot(
##             data=Wealth_cube,
##             aes(
##                 x=time,
##                 y=wealth,
##                 group=interaction(runs,community),
##                 color=community
##             )
##         ) +
##         geom_line(size=0.1) +
##         ylab(expression(log[10]*"(G)")) +
##         ## ylab("log(G)") +
##         xlab("Time") + 
##         ylim(c(-5,5)) + 
##         scale_x_continuous(
##             breaks=seq(0,5000,5000),
##             labels=c("0","50")
##         ) +
##         ## scale_y_continuous(
##         ##     breaks=seq(0,0.5,0.5),
##         ##     limits=c(0,0.5)
##         ## ) +     
##         facet_grid(lambda_P ~ alpha_H, labeller=label_parsed) +
##         scale_color_manual(
##             name=element_blank(),
##             labels=c("Planned","Unplanned"),
##             values=c("cyan","magenta")) + 
##         ## scale_color_discrete() + 
##         theme(
##             aspect.ratio=1,
##             axis.text.x=element_text(size=rel(0.8)),
##             axis.text.y=element_text(size=rel(0.8)),
##             strip.text.x=element_text(size=rel(0.8)),
##             strip.text.y=element_text(size=rel(0.8)),
##             axis.title.x=element_text(size=rel(0.8)),
##             axis.title.y=element_text(size=rel(0.8)),
##             panel.background = element_rect(
##                 fill="transparent",
##                 colour="black",
##                 size=0.25,
##                 linetype=1),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             strip.background=element_blank(),
##             legend.position=c(0.5,-0.1),
##             ## legend.position="bottom",
##             legend.background=element_rect(fill=alpha('white',0)),
##             legend.key=element_rect(fill=NA)
##         ) + guides(color = guide_legend(override.aes = list(alpha = 1), ncol=2, byrow=TRUE))
##         ## theme(plot.margin=add_margin(t=0, unit="cm"))
##     p
## }

## plot_wealth_gini = function(wealth, gini_index, param_combinations, maxtime, dt) {    

##     Gini_cube = 
##         as.tbl_cube(gini_index) %>% 
##         as_tibble %>%
##         mutate_if(sapply(., is.character), as.factor) %>%
##         left_join(param_combinations, by="scenario")    
    
##     Wealth_cube =
##         as.tbl_cube(wealth) %>% 
##         as_tibble %>%
##         mutate_if(sapply(., is.character), as.factor) %>%
##         left_join(param_combinations, by="scenario")

##     data_cube = Gini_cube %>% left_join(Wealth_cube)
    
##     lambda_P_vals =
##         unique(param_combinations$lambda_P) %>%
##         sort %>%
##         format(digits=4, drop0trailing = TRUE)
##     alpha_H_vals =
##         unique(param_combinations$alpha_H) %>%
##         sort %>%
##         format(digits=4, drop0trailing=TRUE)
##     data_cube$lambda_P = factor(
##         data_cube$lambda_P,
##         labels = sapply(lambda_P_vals, FUN=function(x) bquote(lambda[P]==.(x)))
##     )
##     data_cube$alpha_H = factor(
##         data_cube$alpha_H,
##         labels = sapply(alpha_H_vals, FUN=function(x) bquote(alpha[H]==.(x)))
##     )
##     data_cube$community = factor(data_cube$community)
##     data_cube$wealth = log10(data_cube$wealth)
##     dummy_wealth = data_cube %>%
##         dplyr::select(-gini_index) %>%
##         spread(key=community, value=wealth) %>%
##         mutate(`3`=0) %>%
##         gather(community, wealth, `1`, `2`, `3`)
    
##     Gini_cube_limits =
##         data_cube %>%
##         dplyr::select(-wealth) %>% 
##         spread(key=runs, value=gini_index) %>%
##         mutate(
##             gini_max=apply(.[grep("^r", names(.))], 1, max, na.rm=TRUE),
##             gini_min=apply(.[grep("^r", names(.))], 1, min, na.rm=TRUE)
##         ) %>%
##         mutate(
##             gini_max = gini_max * 20 - 5, # scale
##             gini_min = gini_min * 20 - 5  # scale
##         ) %>%
##         mutate(
##             gini_max = ifelse(is.finite(gini_max), gini_max, NA),
##             gini_min = ifelse(is.finite(gini_min), gini_min, NA)
##         )
    
##     Gini_cube_limits$community %<>% as.character

##     ## NEW:
##     Gini_cube_limits %<>% dplyr::select(!starts_with('r'))
    
##     data_cube %<>%
##         left_join(Gini_cube_limits)

##     ## NEW:
##     ## data_cube %<>%
##     ##     right_join(dummy_wealth)

##     p =
##         ## ggplot(
##         ##     data=data_cube,
##         ##     aes(
##         ##         x=time,
##         ##         y=wealth,
##         ##         ymin=gini_min,
##         ##         ymax=gini_max,
##         ##         group=interaction(runs,community),
##         ##         color=community
##         ##     )
##         ## ) +
##         ## geom_ribbon(fill="grey90", alpha=.25, size=0) + #, colour=NA) + 
##         ## geom_line(color=community, size=0.05, alpha=.5) +
##         ggplot(
##             data=data_cube
##         ) +
##         geom_ribbon(
##             aes(ymin=gini_min, ymax=gini_max, x=time, fill="band"), alpha=.25, size=0) +
##         geom_line(aes(x=time, y=wealth, group=interaction(runs,community), color=community), size=0.05, alpha=.75) + 
##         ## geom_ribbon(fill="grey90", alpha=.25, size=0) + #, colour=NA) + 
##         ## geom_line(color=community, size=0.05, alpha=.5) +
##         ylab(expression(log[10]*"(G)")) +
##         xlab("Time") + 
##         scale_x_continuous(
##             breaks=seq(0,maxtime/dt,maxtime/dt),
##             labels=c("0","50")
##         ) +
##         scale_y_continuous(
##             breaks=seq(-5, 5, 10),
##             limits=c(-5,5),
##             sec.axis=sec_axis(
##                 ~(.+5)/20,
##                 breaks=seq(0,0.5,0.5),
##                 name="Gini coefficient"
##             )
##         ) +
##         ## scale_color_manual(
##         ##     name=element_blank(),
##         ##     labels=c("Planned","Unplanned","Gini coefficient"),
##         ##     values=c("magenta","cyan3",NA)
##         ## ) +
##         scale_color_manual(
##             name=element_blank(),
##             labels=c("Planned","Unplanned"),
##             values=c("magenta","cyan3")
##         ) +
##         scale_fill_manual(
##             name=element_blank(),
##             labels=c("Gini coefficient"),
##             values=c("grey50")
##         ) +     
##         ## scale_fill_discrete(guide="none") + 
##         ## facet_grid(mu_S_prime ~ alpha_H, labeller=label_parsed) +
##         facet_grid(lambda_P ~ alpha_H, labeller=label_parsed) +
##         theme(
##             aspect.ratio=1,
##             axis.line = element_line(color="black", size=0.25),
##             axis.text.x=element_text(size=rel(0.8)),
##             axis.text.y=element_text(size=rel(0.8)),
##             strip.text.x=element_text(size=rel(0.8)),
##             strip.text.y=element_text(size=rel(0.8)),
##             axis.title.x=element_text(size=rel(0.8)),
##             axis.title.y=element_text(size=rel(0.8)),
##             panel.background = element_rect(
##                 fill="transparent",
##                 colour="black",
##                 size=0.25,
##                 linetype=1),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             strip.background=element_blank(),
##             legend.position=c(0.5,-0.15),
##             legend.key = element_rect(fill="white"),
##             legend.box = "horizontal"
##         ) +        
##         guides(
##             color=guide_legend(
##                 override.aes=list(
##                     color=c("magenta","cyan3"),
##                     fill="white",
##                     size=0.5,
##                     alpha=c(1,1)),
##                 ncol=3, byrow=TRUE
##             )
##         )
##     p
    
##     ## p =
##     ##     ggplot(
##     ##         data=data_cube,
##     ##         aes(
##     ##             x=time,
##     ##             y=wealth,
##     ##             ymin=gini_min,
##     ##             ymax=gini_max,
##     ##             group=interaction(runs,community),
##     ##             color=community
##     ##         )
##     ##     ) +
##     ##     geom_ribbon(fill="grey90", alpha=.25, size=0) + 
##     ##     geom_line(size=0.05, alpha=.5) +
##     ##     ylab(expression(log[10]*"(G)")) +
##     ##     xlab("Time") + 
##     ##     scale_x_continuous(
##     ##         breaks=seq(0,maxtime/dt,maxtime/dt),
##     ##         labels=c("0","50")
##     ##     ) +
##     ##     scale_y_continuous(
##     ##         breaks=seq(-5, 5, 10),
##     ##         limits=c(-5,5),
##     ##         sec.axis=sec_axis(
##     ##             ~(.+5)/20,
##     ##             breaks=seq(0,0.5,0.5),
##     ##             name="Gini coefficient"
##     ##         )
##     ##     ) +
##     ##     scale_color_manual(
##     ##         name=element_blank(),
##     ##         labels=c("Planned","Unplanned","Gini coefficient"),
##     ##         values=c("magenta","cyan3",NA)
##     ##     ) +
##     ##     scale_fill_discrete(guide="none") + 
##     ##     facet_grid(lambda_P ~ alpha_H, labeller=label_parsed) +
##     ##     theme(
##     ##         aspect.ratio=1,
##     ##         axis.line = element_line(color="black", size=0.25),
##     ##         axis.text.x=element_text(size=rel(0.8)),
##     ##         axis.text.y=element_text(size=rel(0.8)),
##     ##         strip.text.x=element_text(size=rel(0.8)),
##     ##         strip.text.y=element_text(size=rel(0.8)),
##     ##         axis.title.x=element_text(size=rel(0.8)),
##     ##         axis.title.y=element_text(size=rel(0.8)),
##     ##         panel.background = element_rect(
##     ##             fill="transparent",
##     ##             colour="black",
##     ##             size=0.25,
##     ##             linetype=1),
##     ##         panel.grid.major.y=element_blank(),
##     ##         panel.grid.minor.y=element_blank(),
##     ##         strip.background=element_blank(),
##     ##         legend.position=c(0.5,-0.15)
##     ##     ) +        
##     ##     guides(
##     ##         color=guide_legend(
##     ##             override.aes=list(
##     ##                 color=c("magenta","cyan3","grey90"),
##     ##                 fill=c("white","white","grey90"),
##     ##                 alpha=c(1,1,1)), ncol=3, byrow=TRUE
##     ##         )
##     ##     )
##     ## p
## }

## plot_wealth_gini2 = function(wealth, gini_index, param_combinations, maxtime, dt) {    

##     Gini_cube = 
##         as.tbl_cube(gini_index) %>% 
##         as_tibble %>%
##         mutate_if(sapply(., is.character), as.factor) %>%
##         left_join(param_combinations, by="scenario")    
    
##     Wealth_cube =
##         as.tbl_cube(wealth) %>% 
##         as_tibble %>%
##         mutate_if(sapply(., is.character), as.factor) %>%
##         left_join(param_combinations, by="scenario")

##     data_cube = Gini_cube %>% left_join(Wealth_cube)
    
##     mu_S_prime_vals =
##         unique(param_combinations$mu_S_prime) %>%
##         sort %>%
##         format(digits=4, drop0trailing = TRUE)
##     alpha_H_vals =
##         unique(param_combinations$alpha_H) %>%
##         sort %>%
##         format(digits=4, drop0trailing=TRUE)
##     data_cube$mu_S_prime = factor(
##         data_cube$mu_S_prime,
##         labels = sapply(mu_S_prime_vals, FUN=function(x) bquote(paste(plain(mu[S]), plain("'")==.(x))))
##         ## labels = sapply(mu_S_prime_vals, FUN=function(x) bquote(lambda[P]==.(x)))
##     )
##     data_cube$alpha_H = factor(
##         data_cube$alpha_H,
##         labels = sapply(alpha_H_vals, FUN=function(x) bquote(alpha[H]==.(x)))
##     )
##     data_cube$community = factor(data_cube$community)
##     data_cube$wealth = log10(data_cube$wealth)
##     ## dummy_wealth = data_cube %>% dplyr::select(-gini_index) %>% spread(key=community, value=wealth) %>% mutate(`3`=0) %>% gather(community, wealth, `1`, `2`, `3`)
    
##     Gini_cube_limits =
##         data_cube %>%
##         dplyr::select(-wealth) %>% 
##         spread(key=runs, value=gini_index) %>%
##         mutate(
##             gini_max=apply(.[grep("^r", names(.))], 1, max, na.rm=TRUE),
##             gini_min=apply(.[grep("^r", names(.))], 1, min, na.rm=TRUE)
##         ) %>%
##         mutate(
##             gini_max = gini_max * 20 - 5, # scale
##             gini_min = gini_min * 20 - 5  # scale
##         ) %>%
##         mutate(
##             gini_max = ifelse(is.finite(gini_max), gini_max, NA),
##             gini_min = ifelse(is.finite(gini_min), gini_min, NA)
##         )
    
##     Gini_cube_limits$community %<>% as.character
##     Gini_cube_limits %<>% dplyr::select(!starts_with('r'))
    
##     data_cube %<>% left_join(Gini_cube_limits)
##     ## data_cube %<>% right_join(dummy_wealth)
    
##     p =
##         ## ggplot(
##         ##     data=data_cube,
##         ##     aes(
##         ##         x=time,
##         ##         y=wealth,
##         ##         ymin=gini_min,
##         ##         ymax=gini_max,
##         ##         group=interaction(runs,community),
##         ##         color=community
##         ##     )
##         ## ) +
##         ## geom_ribbon(fill="grey90", alpha=.25, size=0) + #, colour=NA) + 
##         ## geom_line(color=community, size=0.05, alpha=.5) +
##         ggplot(
##             data=data_cube
##         ) +
##         geom_ribbon(
##             aes(ymin=gini_min, ymax=gini_max, x=time, fill="band"), alpha=.25, size=0) +
##         geom_line(aes(x=time, y=wealth, group=interaction(runs,community), color=community), size=0.05, alpha=.75) + 
##         ## geom_ribbon(fill="grey90", alpha=.25, size=0) + #, colour=NA) + 
##         ## geom_line(color=community, size=0.05, alpha=.5) +
##         ylab(expression(log[10]*"(G)")) +
##         xlab("Time") + 
##         scale_x_continuous(
##             breaks=seq(0,maxtime/dt,maxtime/dt),
##             labels=c("0","50")
##         ) +
##         scale_y_continuous(
##             breaks=seq(-5, 5, 10),
##             limits=c(-5,5),
##             sec.axis=sec_axis(
##                 ~(.+5)/20,
##                 breaks=seq(0,0.5,0.5),
##                 name="Gini coefficient"
##             )
##         ) +
##         ## scale_color_manual(
##         ##     name=element_blank(),
##         ##     labels=c("Planned","Unplanned","Gini coefficient"),
##         ##     values=c("magenta","cyan3",NA)
##         ## ) +
##         scale_color_manual(
##             name=element_blank(),
##             labels=c("Planned","Unplanned"),
##             values=c("magenta","cyan3")
##         ) +
##         scale_fill_manual(
##             name=element_blank(),
##             labels=c("Gini coefficient"),
##             values=c("grey50")
##         ) +     
##         ## scale_fill_discrete(guide="none") + 
##         facet_grid(mu_S_prime ~ alpha_H, labeller=label_parsed) +
##         theme(
##             aspect.ratio=1,
##             axis.line = element_line(color="black", size=0.25),
##             axis.text.x=element_text(size=rel(0.8)),
##             axis.text.y=element_text(size=rel(0.8)),
##             strip.text.x=element_text(size=rel(0.8)),
##             strip.text.y=element_text(size=rel(0.8)),
##             axis.title.x=element_text(size=rel(0.8)),
##             axis.title.y=element_text(size=rel(0.8)),
##             panel.background = element_rect(
##                 fill="transparent",
##                 colour="black",
##                 size=0.25,
##                 linetype=1),
##             panel.grid.major.y=element_blank(),
##             panel.grid.minor.y=element_blank(),
##             strip.background=element_blank(),
##             legend.position=c(0.5,-0.15),
##             legend.key = element_rect(fill="white"),
##             legend.box = "horizontal"
##         ) +        
##         guides(
##             color=guide_legend(
##                 override.aes=list(
##                     color=c("magenta","cyan3"),
##                     fill="white",
##                     size=0.5,
##                     alpha=c(1,1)),
##                 ncol=3, byrow=TRUE
##             )
##         )
##     p
## }
