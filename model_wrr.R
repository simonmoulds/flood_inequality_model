## Author : Simon Moulds
## Date   : Dec 2019

## Implementation of the socio-hydrological model developed by
## Di Baldassarre et al (2015), published in Water Resources Research

library(ggplot2)

## Set output directory (currently set to current working directory):
outdir = "."

## Time parameters
maxtime = 200
dt = 0.001
t = seq(0, maxtime, dt)
nt = length(t)

## ================================================================= ##
## 1 - Compute inter-arrival time of flood events & flood magnitude  ##
## ================================================================= ##

## In time, this data could be replaced by observed river discharge.

## Because the timing and magnitude of flood events are modelled as a
## random process, here we set the random seed so that each simulation
## is the same (useful for model development). If you remove this line
## you will find that the model produces a different result each time
## (useful for exploring the model space). 
set.seed(42)

t_mean = 1/dt
H_mult = 1
w = rep(0, nt)
index = 0
t0 = 0
t1 = 0
theta=0.28

while (index < nt) {
    x = runif(1)
    t1 = t0 + (t_mean * log(1/(1-x)))
    index = which.min(abs(t/dt - t1))
    y = runif(1)
    w[index] = H_mult * ((theta + 1) / theta) * (1 - (1 - y) ** theta)
    t0 = t1
}

## ================================================================= ##
## 2 - Set model parameters & state variables                        ##
## ================================================================= ##

## See Di Baldassarre et al (2015) for details

## Parameter alpha_H can be used to represent vulnerability (e.g.
## see Ciullo et al. 2017). This could be an entry-point for dynamically
## representing vulnerability, and the impact of trust/vulnerability on
## reducing/increasing.

alpha_H = 10 
xi_H = 0.2
rho_D = 0.03
alpha_D = 5
epsilon_T = 1.1
keta_T = 2e-05
mu_S = 0.06 

## State variables
F = rep(0, nt)
R = rep(0, nt)
D = rep(0.1, nt)
H = rep(0, nt)
M = rep(0, nt)

## ================================================================= ##
## 3 - Run the model                                                 ##
## ================================================================= ##

## N.B. Equation numbers refer to those in Di Baldassarre et al. (2014)

for (i in 1:(nt-1)) {

    ## Equation 1: relative flood damage
    if ((w[i] + xi_H * H[i]) > H[i]) {
        F[i] = 1 - exp(-(w[i] + xi_H * H[i]) / alpha_H)
    } else {
        F[i] = 0
    }
    
    ## Equation 2: increase height of flood protection
    if (F[i] > 0) { 
        R[i] = epsilon_T * (w[i] + xi_H * H[i] - H[i])
    } else {
        R[i] = 0
    }

    ## Value of Dirac comb (always zero except when flooding occurs, when it
    ## assumes a (theoretical) value of positive infinity (i.e. when dt is
    ## infinitely small), with integral equal to one. As we are implementing
    ## a numerical solution we obtain its actual value by dividing by dt
    if (w[i] > 0) {
        Delta = 1/dt
    } else {
        Delta = 0
    }

    ## Equation 3: compute change in state variables
    dD_dt = rho_D * (1 - D[i] * (1 + alpha_D * M[i])) - Delta * F[i] * D[i]
    dH_dt = Delta * R[i] - keta_T * H[i]
    dM_dt = Delta * F[i] * D[i] - mu_S * M[i]

    D[i+1] = D[i] + dD_dt * dt
    H[i+1] = H[i] + dH_dt * dt
    M[i+1] = M[i] + dM_dt * dt
}

## ================================================================= ##
## 4 - Plot model results                                            ##
## ================================================================= ##

## Create a time series data frame with model state variables
output = data.frame(
    time=t,
    density=D,
    level=H,
    memory=M,
    damage=F
)

## Set some aspects of the plot theme (see ggplot2 documentation for
## more information about this)
theme_set(
    theme_bw() +
    theme(
        axis.text.x=element_text(size=rel(0.8)), # size of text used for x-axis label
        axis.text.y=element_text(size=rel(0.8)), # size of text used for y-axis label
        ## strip.background=element_blank(),        
        ## strip.text.x=element_blank(),
        axis.title.x=element_text(size=rel(0.8)), # size of text used for x-axis title
        axis.title.y=element_text(size=rel(0.8)), # size of text used for y-axis title
        panel.grid.major.y=element_blank(),       # remove major/minor x/y grid lines
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank()
    )
)

## Plot population density
p1 =
    ggplot(
        data=output,
        aes(x=time, y=density)
    ) +   
    geom_line(size=0.3, alpha=1) +
    ylab("D") +
    xlab(element_blank()) +
    ylim(c(0,1)) +
    scale_x_continuous(
        breaks=seq(0,50/dt,10/dt),
        labels=c("0","10","20","30","40","50")
    ) +
    labs(tag="(a)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        strip.text.x=element_text(size=rel(0.8)),
        legend.position="none"
    )
p1

## Plot flood protection level
p2 =
    ggplot(
        data=output,
        aes(x=time, y=level)
    ) +                    
    geom_line(size=0.3, alpha=1) +
    ylab("H") +
    xlab(element_blank()) +
    ylim(c(0,5)) + 
    scale_x_continuous(
        breaks=seq(0,50/dt,10/dt),
        labels=c("0","10","20","30","40","50")
    ) +
    labs(tag="(b)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="none"
    )

## Plot relative flood damage
p3 =
    ggplot(
        data=output,
        aes(x=time, y=damage)
    ) +    
    geom_point(size=0.75, alpha=1) +
    geom_line(size=0.3, alpha=0) + 
    ylab("F") +
    ylim(c(0,1)) + 
    scale_x_continuous(
        breaks=seq(0,50/dt,10/dt),
        labels=c("0","10","20","30","40","50")
    ) +
    labs(tag="(c)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="bottom"
    )

## Plot memory of flood events
p4 =
    ggplot(
        data=output,
        aes(x=time, y=memory)
    ) +    
    geom_point(size=0.75, alpha=1) +
    geom_line(size=0.3, alpha=0) + 
    ylab("F") +
    xlab("Time") +
    ylim(c(0,1)) + 
    scale_x_continuous(
        breaks=seq(0,50/dt,10/dt),
        labels=c("0","10","20","30","40","50")
    ) +
    labs(tag="(d)") + 
    theme(
        plot.tag=element_text(size=rel(0.8)),
        legend.position="bottom"
    )

## Join p1-4
pl <- lapply(list(p1,p2,p3,p4), ggplotGrob)
pl <- Reduce(gridExtra::gtable_rbind, pl)

## Display the combined plots
grid::grid.draw(pl)

## Write output to file
png(
    file.path(outdir, "wrr_state_vars.png"),
    width=6,
    height=8,
    units="in",
    res=320
)
grid::grid.draw(pl)
dev.off()
