# Area level Fay Herriot 

The classical area-level Fay_herriot model is defined by:

$$\hat{\theta}_d | \theta_{d} \sim N(\theta_{d}, \psi_d^2)$$
The above say that our district survey estimates are a noisy measurement of the true underlying district mean $\theta_{d}$, where $\psi_{d}^{2}$ is the known sampling variance from the survey design (measurement error stage)

$$\theta_{d} = X_{d}^{T}\beta + \mu_{d}$$
The true district mean is then a function of district level covarites $x_{d}$ plus a district random effect $u_{d}$. Everything is already aggregated to the district, there are no clusters, no children, no individual observations.

$$u_{d} \stackrel{\text{iid}}{\sim} N(0, \sigma_{u}^{2})$$
The random effects are independent across districts, there is no spatial structure, no borrowing from neighbors, just shrinkage towards the co variate-predicated mean



# The Multilevel Fay Heroit Model

##Sampling Model

The sampling model is similar for all models to be defined iin the subsequent sections as they all share the first layer. The smapling model describes how the data arise from the underlying true malaria risk

$$ y_{c} | n_{c}, \pi_{c} \sim Binomial(n_{c}, \pi_{c}) $$
$$logit(\pi_{c}) = \eta_{c}$$

The sampling model above states that each child tested is either positive or negative. With $n_{c}$ children tested in cluster c, the
number of positives follows a Binomial distribution governed by the true underlying prevalence $\pi_{c}$. We model the log-odds of prevalence
rather than the probability directly, because log-odds are unbounded and work naturally with linear models. And the quantity $\eta_{c}$ is what each model below specifies

# Model 1 - Multilevel IID (With no covariates)

The IID model with no covariates ask if there is meaningful variation in malaria risk across districts after accounting for the fact that clusters within the same district tend to look the same (correlated). The model, which in this
case is the linking model (layer of the fay-herroit) is defined as:

$$\eta_{c} = \mu + \sigma_{u}u_{d[c]} + \sigma_{v}\epsilon_{c}$$

Where: \n
$\mu$ is the national baseline. The garage log-odds of malaria across all of Ghana. Every cluster starts here. That is if we randomly picked a child from a random cluster in Ghana, what would we expect, before we account for district or cluster differences.

$u_{d[c]}$ is the district random effect. A positive or negative adjustment for the district that each cluster c is found in. Every cluster in the same district share this same value, have this single number. A district with high malaria burden gets a positive $u_{d}$ and the opposite is true
. A because the prior we are going to give to this random effect is centered around zero, districts with little or no data get pulled back toward the national average and this is the borrowing of strength we are after in SAE. 

$\sigma_{u}$ between district variation.  This defines the scaling weight on the district effect. If it is large it means districts in Ghana qenuinely differ a great deal from each other. If it is approxiimately zero, the all districts
are essentially the same and the district effect contributes nothing. This is estimated from the data directly

$\epsilon_{c}$ cluster deviation. Even within the same district, individual clusters differ from each other. One cluster in Brong-Ahafo might be near a stagnant lake; another in the same district sits on higher ground. This cluster-to-cluster wobble within a district is $\epsilon_{c}$
$\epsilon_{c}$. It captures local heterogeneity that the district-level effect cannot see.

$\sigma_{v}$ within-district cluster variation. The scaling weight on the cluster deviation. This tells you how noisy and heterogeneous clusters are even within the same district. Statistically, this absorbs overdispersion, the extra variation that exceeds what a simple Binomial model expects. Without this term, uncertainty intervals would be overconfident.

Total variation from the multilevel IID is the explained by:

$$Total Variation in \eta_{c} = \sigma_{u}^{2} + \sigma_{v}^{2}$$

If $\sigma_{u}^{2}$ dominates: malaria is mostly a district-level story. Knowing which district you are in explains most of the risk. SAE is very effective.

If $\sigma_{v}^{2}$ dominates: districts are internally heterogeneous. Even knowing the district does not help much for a specific cluster. SAE is less powerful but uncertainty is correctly represented.