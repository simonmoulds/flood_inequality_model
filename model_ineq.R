## Author : Simon Moulds
## Date   : Dec 2019 - Oct 2020

library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(egg)
library(cubelyr)

source("model_funs.R")

## CHANGE THIS IF RUNNING ON A DIFFERENT MACHINE
outdir = "~/Dropbox/plots"

## Parameter description:
## ----------------------
## Symbol    | Description
## ---------------------------------------------------------
## xi_H      | proportion of additional high water level due
##           |  to levee heightening
## gamma_E   | cost for unit height R and width sqrt(G) of
##           | levee raising (0.005 -> inf)
## epsilon_T | safety factor for raising of levees
## keta_T    | rate of decay of levees
## alpha_H   | slope of floodplain, resilience of human
##           | settlement
## rho_E     | maximum relative growth rate 
## lambda_P  | distance at which people would accept to live,
##           | remembering past floods (0.005 -> inf)
## phi_P     | rate by which new properties can be built
## alpha_S   | proportion of shock after flooding if levees
##           | are raised (0 -> 1)
## mu_S      | memory loss rate (0.01 -> 10)
## tau       | redistribution

## Initial values of state variables

## name      | D-B et al. (2013) | Viglione et al. (2014)
## ------------------------------------------------------
## F_0       | 0                 | 0
## G_0       | 100               | 0.01
## D_0       | 2000              | 1
## R_0       | 0                 | 0
## H_0       | 0                 | 0
## S_0       | 0                 | 0
## M_0       | 0                 | 0

## Parameter values

## name      | D-B et al. (2013) | Viglione et al. (2014) | Blair (2017) [MRes thesis]
## ---------------------------------------------------------------------
## xi_H      | 0.5               | 0.5                    | 0.2
## gamma_E   | 0.5, 5, 5000      | 0.005                  | 
## epsilon_T | 1.1               | 1.1                    | 1.1
## keta_T    | 0.003 0.0003      | 0.1                    | 0.0002
## alpha_H   | 0.1               | 10                     | 10
## rho_E     | 0.02              | 1                      | 
## lambda_P  | 12000             | 0-5                    | 
## phi_P     | 100               | 0.1                    | 
## alpha_S   | 0.5               | 0.5                    | 
## mu_S      | 0.05              | 0.01-10                | 0.12
## lambda_E  | 5000              | -                      |
## keta_W    |                   |                        |

## New parameters:
## ===============
## tau       | 0, 0.5

## ################################## ##
## ################################## ##
##                                    ##
## Experiment 1: Sensitivity analysis ##
##                                    ##
## ################################## ##
## ################################## ##

## Set up the simulation
## #####################

## Reduce dt, otherwise the size of the resulting data frame
## becomes too large and R will crash
maxtime = 50
dt = 0.1
H_mult = 1
theta = 0.28

## Number of communities/scenarios/time steps/model runs
nc = 2
ns = 36
nt = maxtime / dt + 1
nr = 50
## nr = 1

## Scenario 1: Build back better, leave the poor behind
## ####################################################

## "High efficiency, low equity"

## Key points:

## * Planned settlement installs flood protection, while unplanned settlement cannot
## * Economic imperative means that unplanned settlement continues to move
##   towards places which enable economic activity
## * Prosperity in planned settlement does "trickle down" to an extent

## High efficiency is enabled by setting the safety factor associated with flood
## protection to 1.1 for the planned settlement, enabling this community to
## install adequate flood protection following a flood event. 

## initial values of state variables
init = get_default_initial_values()

## default parameter set
params = get_default_parameter_set(ns, nc)

params$mu_S[] = 0.25
params$rho_E[] = 1                           # economic growth
params$epsilon_T[] = rep(c(1.1, 0), each=ns) # levee safety factor
params$tau[] = 0                             # wealth redistribution

## We then select two parameters to vary (lambda_P, alpha_H)
param_grid = expand.grid(
    lambda_P = seq(1/3, 2, length.out=6),
    alpha_H = seq(4, 10, length.out=6)
)
params$lambda_P[,1] = 2
params$lambda_P[,2] = param_grid$lambda_P
params$alpha_H[,1] = 10
params$alpha_H[,2] = param_grid$alpha_H

## label each combination of lambda_P/alpha_H for
## the unplanned settlement, to help with plots later
param_combinations = data.frame(
    scenario=paste0("s", 1:36),
    lambda_P = param_grid$lambda_P,
    alpha_H = param_grid$alpha_H
)

## Create arrays to hold results
gini_index = array(data=NA, dim=c(ns, nt, nr))
dimnames(gini_index) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "runs" = sprintf("r%d", 1:nr)
)
wealth = array(data=NA, dim=c(ns, nt, nc, nr))
dimnames(wealth) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "community" = 1:nc,
    "runs" = sprintf("r%d", 1:nr)
)

## Run simulations (set random seed to ensure reproducibility)
set.seed(42)
for (i in 1:nr) {
    print(i)
    w = get_water_level_ts(maxtime=maxtime, dt=dt)
    out = model_ineq(w, init, params, nc, ns)
    wealth[,,,i] = out$G
    gini_index[,,i] = compute_gini(out$G)
}

## Plot results
p1 = plot_wealth_gini(wealth, gini_index, param_combinations, maxtime, dt)
ggsave(
    file.path(outdir, "flood_ineq_gini_scen1.png"),
    p1, width=6, height=7, units="in"
)

## Scenario 2: Build back worse, leave the poor behind
## ###################################################

## "Low efficiency, low equity"

## We represent low efficiency as low economic growth
params$mu_S[] = 0.25
params$rho_E[] = 0.5
params$epsilon_T[] = rep(c(1.1, 0), each=ns)
params$tau[] = 0

## Preallocate arrays to store results
gini_index = array(data=NA, dim=c(ns, nt, nr))
dimnames(gini_index) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "runs" = sprintf("r%d", 1:nr)
)

wealth = array(data=NA, dim=c(ns, nt, nc, nr))
dimnames(wealth) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "community" = 1:nc,
    "runs" = sprintf("r%d", 1:nr)
)

## Run simulations:
set.seed(42)
for (i in 1:nr) {
    print(i)
    w = get_water_level_ts(maxtime=maxtime, dt=dt)
    out = model_ineq(w, init, params, nc, ns)
    wealth[,,,i] = out$G
    gini_index[,,i] = compute_gini(out$G)
}

## Make plots
p1 = plot_wealth_gini(wealth, gini_index, param_combinations, maxtime, dt)

ggsave(
    file.path(outdir, "flood_ineq_gini_scen2.png"),
    p1, width=6, height=7, units="in"
)

## Scenario 3: Build back worse, leave no one behind
## #################################################

## "low efficiency, high equity"

## Preallocate arrays to store results:
params$mu_S[] = 0.25
params$rho_E[] = 0.5
params$epsilon_T[] = rep(c(1.1, 1.1), each=ns)
params$tau[] = 0.5
params$mu_S_prime[] = 2

gini_index = array(data=NA, dim=c(ns, nt, nr))
dimnames(gini_index) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "runs" = sprintf("r%d", 1:nr)
)

wealth = array(data=NA, dim=c(ns, nt, nc, nr))
dimnames(wealth) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "community" = 1:nc,
    "runs" = sprintf("r%d", 1:nr)
)

## Run simulations:
set.seed(42)
for (i in 1:nr) {
    print(i)
    w = get_water_level_ts(maxtime=maxtime, dt=dt)
    out = model_ineq(w, init, params, nc, ns)
    wealth[,,,i] = out$G
    gini_index[,,i] = compute_gini(out$G)
}

## Plot results
p1 = plot_wealth_gini(wealth, gini_index, param_combinations, maxtime, dt)
ggsave(
    file.path(outdir, "flood_ineq_gini_scen3.png"),
    p1, width=6, height=7, units="in"
)

## Scenario 4: Build back better, leave no one behind
## ##################################################

## "high efficiency, high equity"

params$mu_S[] = 0.25
params$rho_E[] = 1
params$epsilon_T[] = rep(c(1.1, 1.1), each=ns)
params$mu_S_prime[] = 2
params$tau[] = 0.5

gini_index = array(data=NA, dim=c(ns, nt, nr))
dimnames(gini_index) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "runs" = sprintf("r%d", 1:nr)
)

wealth = array(data=NA, dim=c(ns, nt, nc, nr))
dimnames(wealth) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "community" = 1:nc,
    "runs" = sprintf("r%d", 1:nr)
)

## Run simulations:
set.seed(42)
for (i in 1:nr) {
    print(i)
    w = get_water_level_ts(maxtime=maxtime, dt=dt)
    out = model_ineq(w, init, params, nc, ns)
    wealth[,,,i] = out$G
    gini_index[,,i] = compute_gini(out$G)
}

## Make plots
p1 = plot_wealth_gini(wealth, gini_index, param_combinations, maxtime, dt)
ggsave(
    file.path(outdir, "flood_ineq_gini_scen4.png"),
    p1, width=6, height=7, units="in"
)

## ################################## ##
## ################################## ##
##                                    ##
## Experiment 2: System archetypes    ##
##                                    ##
## ################################## ##
## ################################## ##

maxtime = 50
dt = 0.01
H_mult = 1
theta = 0.28

nc = 2                # number of communities (planned, unplanned)
ns = 4                # number of scenarios (= number of system archetypes)
nt = maxtime / dt + 1 # number of time steps
nr = 1                # number of model runs

## initial values of state variables
init = get_default_initial_values()
## default parameter set
params = get_default_parameter_set(ns, nc)

## Scenario 1: Build back better, leave the poor behind
## ####################################################

## "High efficiency, low equity"

params$rho_E[1,] = c(1, 1)
params$epsilon_T[1,] = c(1.1, 0)
params$lambda_P[1,] = c(2, 0.666666)
params$alpha_H[1,] = c(10, 5.2)
params$tau[1,] = c(0.00, 0.00)
params$mu_S_prime[] = 10

## Scenario 2: Build back worse, leave the poor behind
## ###################################################

## "Low efficiency, low equity"

params$rho_E[2,] = c(0.5,0.5)
params$epsilon_T[2,] = c(1.1, 0)
params$lambda_P[2,] = c(2, 0.666666)
params$alpha_H[2,] = c(10, 5.2)
params$tau[2,] = c(0, 0)
params$mu_S_prime[] = 10

## Scenario 3: Build back worse, leave no one behind
## #################################################

## "low efficiency, high equity"

params$rho_E[3,] = c(0.5, 0.5)
params$epsilon_T[3,] = c(1.1, 1.1)
params$lambda_P[3,] = c(2, 0.666666)
params$alpha_H[3,] = c(10, 5.2)
params$tau[3,] = c(0.5, 0.5)
params$mu_S_prime[] = 2

## Scenario 4: Build back better, leave no one behind
## ##################################################

## "high efficiency, high equity"

params$rho_E[4,] = c(1, 1)
params$epsilon_T[4,] = c(1.1, 1.1)
params$lambda_P[4,] = c(2, 0.666666)
params$alpha_H[4,] = c(10, 5.2)
params$tau[4,] = c(0.5, 0.5)
params$mu_S_prime[] = 2

## Create arrays to hold results
gini_index = array(data=NA, dim=c(ns, nt, nr))
dimnames(gini_index) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "runs" = sprintf("r%d", 1:nr)
)
wealth = array(data=NA, dim=c(ns, nt, nc, nr))
dimnames(wealth) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "community" = 1:nc,
    "runs" = sprintf("r%d", 1:nr)
)
level = wealth
distance = wealth
damage = wealth

## Run model
set.seed(42)
for (i in 1:nr) {
    print(i)
    w = get_water_level_ts(maxtime=maxtime, dt=dt)
    out = model_ineq(w, init, params, nc, ns)
    wealth[,,,i] = out$G
    level[,,,i] = out$H
    distance[,,,i] = out$D
    damage[,,,i] = out$F
    gini_index[,,i] = compute_gini(out$G)
}

labels=c(
    "High efficiency, low equity ",
    "Low efficiency, low equity  ",
    "Low efficiency, high equity ",
    "High efficiency, high equity"
)

## Make plot
G_cube =
    as.tbl_cube(wealth) %>%
    as_tibble %>%
    mutate_if(sapply(., is.character), as.factor) %>%
    mutate_at(vars(matches("community")), as.factor)
levels(G_cube$community) = c("Planned","Unplanned")

H_cube =
    as.tbl_cube(level) %>%
    as_tibble %>%
    mutate_if(sapply(., is.character), as.factor)
levels(H_cube$community) = c("Planned","Unplanned")

D_cube =
    as.tbl_cube(distance) %>%
    as_tibble %>%
    mutate_if(sapply(., is.character), as.factor)
levels(D_cube$community) = c("Planned","Unplanned")

F_cube =
    as.tbl_cube(damage) %>%
    as_tibble %>%
    mutate_if(sapply(., is.character), as.factor)
levels(F_cube$community) = c("Planned","Unplanned")
F_cube[F_cube == 0] = NA

theme_set(
    theme_bw() +
    theme(
        axis.text.x=element_text(size=rel(0.8)),
        axis.text.y=element_text(size=rel(0.8)),
        strip.background=element_blank(),
        strip.text.x=element_blank(),
        axis.title.x=element_text(size=rel(0.8)),
        axis.title.y=element_text(size=rel(0.8)),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank()
    )
)

p1 =
    ggplot(
        data=G_cube,
        aes(x=time,
            y=log10(wealth),
            group=interaction(runs,scenario),
            color=scenario
            )
    ) +   
    geom_line(size=0.3, alpha=1) +
    ylab(expression(log[10]*"(G)")) +
    xlab(element_blank()) +
    ylim(c(-5,5)) +
    scale_x_continuous(breaks=seq(0,5000,1000), labels=c("0","10","20","30","40","50")) +
    facet_wrap(.~community) + 
    labs(tag="(a)") + 
    theme(
        ## plot.tag.position=c(0.08,0.99),
        plot.tag=element_text(size=rel(0.8)),
        strip.text.x=element_text(size=rel(0.8)),
        legend.position="none"
    )
p1

p2 =
    ggplot(
        data=H_cube,
        aes(x=time,
            y=level,
            group=interaction(runs,scenario),
            color=scenario
            )
    ) +                    
    geom_line(size=0.3, alpha=1) +
    ylab("H") +
    xlab(element_blank()) +
    ylim(c(0,5)) + 
    scale_x_continuous(
        breaks=seq(0,5000,1000),
        labels=c("0","10","20","30","40","50")
    ) +
    facet_wrap(.~community) +
    labs(tag="(b)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="none"
    )

p3 =
    ggplot(
        data=D_cube,
        aes(x=time,
            y=distance,
            group=interaction(runs,scenario),
            color=scenario
            )
    ) +
    geom_line(size=0.3, alpha=1) +
    ylab("D") +
    xlab(element_blank()) +
    ylim(c(0,2.5)) + 
    scale_x_continuous(
        breaks=seq(0,5000,1000),
        labels=c("0","10","20","30","40","50")
    ) +
    facet_wrap(.~community) + 
    labs(tag="(c)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="none"
    )

levels(F_cube$scenario) = labels
p4 =
    ggplot(
        data=F_cube,
        aes(x=time,
            y=damage,
            group=interaction(runs,scenario),
            color=scenario
            )
    ) +    
    geom_point(size=0.75, alpha=1) +
    geom_line(size=0.3, alpha=0) + 
    ylab("F") +
    xlab("Time") +
    ylim(c(0,1)) + 
    scale_x_continuous(breaks=seq(0,5000,1000), labels=c("0","10","20","30","40","50")) +
    facet_wrap(.~community) + 
    scale_color_discrete(name=element_blank()) + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="bottom"
    ) +
    labs(tag="(d)") + 
    guides(
        color=guide_legend(
            override.aes=list(
                alpha=1,
                fill=c("white","white","white","white")
            ),
            ncol=2,
            byrow=TRUE
        )
    ) 

pl <- lapply(list(p1,p2,p3,p4), ggplotGrob)
pl <- Reduce(gridExtra::gtable_rbind, pl)

grid::grid.draw(pl)

png(
    file.path(outdir, "flood_ineq_state_vars.png"),
    width=6,
    height=8,
    units="in",
    res=320
)
grid::grid.draw(pl)
dev.off()

## Plot Gini coefficient and combined wealth
Gini_cube =
    as.tbl_cube(gini_index) %>% 
    as_tibble %>%
    mutate_if(sapply(., is.character), as.factor)

p1 =
    ggplot(
        data=Gini_cube,
        aes(x=time,
            y=gini_index,
            group=interaction(runs,scenario),
            color=scenario)
    ) +
    geom_line(size=0.3, alpha=1) +
    ylab("Gini coefficient") +
    xlab("Time") +
    ylim(c(0,0.5)) +
    scale_x_continuous(
        breaks=seq(0,5000,1000),
        labels=c("0","10","20","30","40","50")
    ) +
    scale_color_discrete(
        name=element_blank(),
        labels=labels
    ) +
    labs(tag="(a)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="bottom"
    ) +
    guides(
        color = guide_legend(
            override.aes = list(alpha = 1),
            ncol=2,
            byrow=TRUE
        )
    )

## Add up wealth from two communities
wealth_combined = apply(wealth, c(1,2,4), FUN=sum, na.rm=TRUE)
wealth_cube =
    as.tbl_cube(wealth_combined) %>% 
    as_tibble %>%
    mutate_if(sapply(., is.character), as.factor) %>%
    mutate_at("wealth_combined", list(wealth_combined=log))

p2 =
    ggplot(
        data=wealth_cube,
        aes(x=time,
            y=wealth_combined,
            group=interaction(runs,scenario),
            color=scenario)
    ) +
    geom_line(size=0.3, alpha=1) +
    ylab(expression(log[10]*"(G)")) +
    xlab("Time") +
    scale_x_continuous(
        breaks=seq(0,5000,1000),
        labels=c("0","10","20","30","40","50")
    ) +
    scale_color_discrete(
        name=element_blank(),
        labels=labels
    ) +
    labs(tag="(b)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="bottom"
    ) +
    guides(
        color = guide_legend(
            override.aes = list(alpha = 1),
            ncol=2,
            byrow=TRUE
        )
    )

legend_b <- get_legend(
    p1 + theme(legend.direction = "horizontal",
               legend.justification="center",
               legend.box.just = "bottom"
               )
)

prow = plot_grid(
    p1 + theme(aspect=1, legend.position="none"),
    p2 + theme(aspect=1, legend.position="none")
)

p = plot_grid(prow, legend_b, ncol=1, rel_heights=c(1,0.2))

png(
    file.path(outdir, "flood_ineq_gini_coef_wealth.png"),
    width=6,
    height=4,
    units="in",
    res=320
)
grid::grid.draw(p)
dev.off()

## ################################## ##
## ################################## ##
##                                    ##
## Experiment 3: Society memory       ##
##                                    ##
## ################################## ##
## ################################## ##

maxtime = 50
dt = 0.1
H_mult = 1
theta = 0.28

nc = 2                # number of communities (planned, unplanned)
ns = 36               # number of scenarios
nt = maxtime / dt + 1 # number of time steps
nr = 50               # number of model runs

## initial values of state variables
init = get_default_initial_values()
## default parameter set
params = get_default_parameter_set(ns, nc)

params$lambda_P[] = rep(c(2, 0.666666), each=ns)
params$mu_S_dash[] = 2
params$tau[] = 0.5

## We then select two parameters to vary (lambda_P, alpha_H)
param_grid = expand.grid(
    mu_S_prime = seq(0, 10, length.out=6),
    alpha_H = seq(4, 10, length.out=6)
)
params$mu_S_prime[,1] = param_grid$mu_S_prime
params$mu_S_prime[,2] = param_grid$mu_S_prime
params$alpha_H[,1] = 10
params$alpha_H[,2] = param_grid$alpha_H

## label each combination of lambda_P/alpha_H for
## the unplanned settlement, to help with plots later
param_combinations = data.frame(
    scenario=paste0("s", 1:36),
    mu_S_prime = param_grid$mu_S_prime,
    alpha_H = param_grid$alpha_H
)

## Create arrays to hold results
gini_index = array(data=NA, dim=c(ns, nt, nr))
dimnames(gini_index) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "runs" = sprintf("r%d", 1:nr)
)
wealth = array(data=NA, dim=c(ns, nt, nc, nr))
dimnames(wealth) = list(
    "scenario" = sprintf("s%d", 1:ns),
    "time" = 1:nt,
    "community" = 1:nc,
    "runs" = sprintf("r%d", 1:nr)
)

set.seed(42)
for (i in 1:nr) {
    print(i)
    w = get_water_level_ts(maxtime=maxtime, dt=dt)
    out = model_ineq(w, init, params, nc, ns)
    wealth[,,,i] = out$G
    gini_index[,,i] = compute_gini(out$G)
}

## Plot results
source("model_funs.R")
p1 = plot_wealth_gini2(wealth, gini_index, param_combinations, maxtime, dt)
ggsave(
    file.path(outdir, "flood_ineq_gini_scen4_mem.png"),
    p1, width=6, height=7, units="in"
)

