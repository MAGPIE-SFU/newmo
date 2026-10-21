# Extrapolating carriage prevalence from invasiveness data

This document details the method used to extrapolate carriage prevalence from invasive prevalence and estimated invasiveness parameters from Nick's paper.

## Overview

There is information about the carriage prevalence from two sources:

- indirectly from invasivness data, which we have more of, and
- directly from carriage data, which we have less of.

We are estimating $x_{i,g}(t)$, the carriage prevalence of serotype $i$ in group $g$ for year $t$ in the child population.

The outline of the procedure is:

1. Extrapolate counts of invasive cases to carriage prevalence using the invasiveness quantities estimated in Nick's paper

1. Combine information from child, adult, and unknown populations by taking a weighted average of the three

1. Define a Dirichlet-multinomial posterior distribution with:

  - observed carriage counts as the mutinomial data and

  - extrapolated prevalence (from the invasive data) and a concentration hyperparameter defining the Dirichlet prior.


## Extrapolating carriage prevalence from invasive data

These are the sources for the extrapolation

| Values | Definition | Source |
| --- | --- | --- |
| $v_s$ | Invasiveness | Nick's paper |
| $I_{s,g}(t)$ | No. invasive samples | SA data |

We extrapolate the number of invasive samples to the number of carriage samples by dividing by Nick's invasiveness:
$$C^C_{s,g}(t) = I_{s,g}(t) / v_s,$$
where the $C$ in the exponent stands for the child population. Similar quantities are defined for the adult ($A$) and unknown ($U$) populations.

This gives a corresponding prevalence estimate for each serotype/group pairing in the child population:
$$p^C_{s,g}(t) = \sum_{s,g}C_{s,g}(t),$$
and similarly for the other two populations.


Correspondence with between variables here and variables in the code:
| In-text variable | In-code variable |
| --- | --- |
| $I_{s,g}(t)$ | `ipd_n` |
| $C^C_{s,g}(t)$ | `r_C` |
| $C^A_{s,g}(t)$ | `r_A` |
| $C^U_{s,g}(t)$ | `r_U` |
| $p^C_{s,g}(t)$ | `p_C_inv` |
| $p^A_{s,g}(t)$ | `p_A_inv` |
| $p^U_{s,g}(t)$ | `p_U_proxy` |



## Constructing the Dirichlet prior

The Dirichlet prior has two parameters: a mean vector of proportions and a concentration parameter.

The prior mean is defined as a weighted average of the extrapolated prevalences, based on fixed weights defined below
$$\tilde p_{s,g}(t) = (1-\tau_A-\tau_U)p_{s,g}^C(t) + \tau_Ap_{s,g}^A(t) + \tau_Up_{s,g}^U(t)$$

The prior parameters are fixed beforehand:
| Hyperparameter | Definition | Current value |
| --- | --- | --- |
| $k$ | Dirichlet concentration | 100 |
| $\tau_A$ | Weight on adult prevalence | 0 |
| $\tau_U$ | Weight on unknown prevalence | 0.0518 |


### Adjusting for zero counts

Since some of the serotype/group combinations have counts of 0 for some years, the formula above needs to be modified to account for this.
This modification just explicitly sets $\tau_A$ or $\tau_U$ to 0 when those respective populations are absent for a particular serotype/group/year combination.
The rest of the algebra in this subsection is just adjusting the weights so that the weighted average is still reasonable.

First, re-define the adult and unknown weights to explicitly account for zeros:
```math
\tau_A^*=\begin{cases} \tau_A\text{,~~}\sum_{s,g}p_{s,g}^A(t)>0 \\ 0\text{,~~~~otherwise}\end{cases}
```
and similarly for $\tau_U^*$. We then take $\tau_{\text{total}}=1-(\tau_A-\tau_A^*) - (\tau_U-\tau_U^*)$.
And then the weights are rescaled:
```math
\begin{aligned}
w_C &= (1-\tau_A-\tau_U)/\tau_\text{total}\\
w_A &= \tau_A^*/\tau_\text{total}\\
w_U &= \tau_U^*/\tau_\text{total}
\end{aligned}
```

So the zero-adjusted prior mean is:
$$\tilde p_{s,g}(t) = w_Cp_{s,g}^C(t) + w_Ap_{s,g}^A(t) + w_Up_{s,g}^U(t)$$


## Posterior estimates

This posterior distribution is Dirichlet, so it has a closed form for the posterior mean.
From the observed carriage data we have the number of carriage samples for each serotype/group/year combination, $n_{s,g}(t)$, and the total number of carriage samples in a particular year, $N_(t)$. This gives us the fraction $\omega(t)=\frac{N(t)}{N(t)+k}$, which can be interpreted as a kind of evidence weight.
The posterior mean is then:
$$\hat p_{s,g}(t) = \omega(t)\frac{n_{s,g}(t)}{N(t)}+(1-\omega(t))\tilde p_{s,g}(t).$$

And the current estimate of the standard deviation of $\hat p_{s,g}(t)$ is

$$SD(\hat p_{s,g}(t)) = \sqrt{\frac{\hat p_{s,g}(t)(1-\hat p_{s,g}(t))}{(N(t)+k+1)}}$$


The lookup table for variables here and in the code is:

| In-text variable | In-code variable |
| --- | --- |
| $\tau_A^*$ | `tauA_eff` |
| $\tau_U^*$ | `tauU_eff` |
| $\tilde p_{s,g}(t)$ | `p_C_aug` |
| $n_{s,g}(t)$ | `n_carriage` |
| $N(t)$ | `N_car` |
| $\hat p_{s,g}(t)$ | `p_final` |
| $SD(\hat p_{s,g}(t))$ | `p_final_sd` |
