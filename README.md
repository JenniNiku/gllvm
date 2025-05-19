Let the marginal log-likelihood for a logistic GLLVM be:

with
*η*<sub>*i**j*</sub> = *β*<sub>0*j*</sub> + **z**<sub>*i*</sub><sup>⊤</sup>**γ**<sub>*j*</sub>.
The above expectation of the first term does not possess a closed-form
solution. [Akin to Gibbs sampling of logistic
regression](https://arxiv.org/pdf/1205.0310), we apply a polya-gamma
data augmentation and additional introduce the auxiliary variable
$\omega \sim PG(b,c) = \frac{\exp(-c^2\omega/2)f(\omega;b,0)}{\mathbb{E}\\\exp(-\frac{1}{2}c^2\omega)\\}$,
[a polya-gamma random
variable](https://stephens999.github.io/fiveMinuteStats/polya_gamma.html),
so that the joint log-likelihood instead becomes:

with corresponding marginal log-likelihood:

Applying Jensen’s inequality, we have the typical VA ELBO:

where
$\mathbb{E}\_{\textbf{z}\_i}\\\log p(y\_{ij}\vert \omega\_{ij}, \textbf{z}\_i)\\ = -\log{2} + (y\_{ij} - \frac{1}{2})\mathbb{E}\_{\textbf{z}\_i}(\eta\_{ij}) -\frac{1}{2}\mathbb{E}\_{\omega\_{ij}}(\omega)\mathbb{E}\_{\textbf{z}\_i}(\eta\_{ij}^2)$
under the assumption
*q*(*ω*<sub>*i**j*</sub>, **z**<sub>*i*</sub>) = *q*(*ω*<sub>*i**j*</sub>)*q*(**z**<sub>*i*</sub>).

<!-- Note, that we cannot proceed under the assumption $q(\omega_{ij},\textbf{z}_i) = q(\omega_{ij}|\textbf{z}_i)q(\textbf{z}_i)$, as it would leave us with $\mathbb{E}_{\textbf{z}_i}\{\log\cosh(\vert\eta\vert/2)\}$ in the likelihood, which does not possess a closed-form solution. -->

it is straightforward to show that if we formulate
*q*(**z**<sub>*i*</sub>|*ω*<sub>*i**j*</sub>)*q*(*ω*<sub>*i**j*</sub>),
the first distribution is conditionally normal and thus possesses a
closed form solution w.r.t. integration of *ω*<sub>*i**j*</sub>. Thus,
due to the polya-gamma augmentation, so we can assume (as usual)
*q*(**z**<sub>*i*</sub>) = 𝒩(**a**<sub>*i*</sub>, **A**<sub>*i*</sub>)
where **a**<sub>*i*</sub> and **A**<sub>*i*</sub> are free variational
parameters to estimate. For *q*(*ω*<sub>*i**j*</sub>):

which we recognize as the (unnormalized) exponentially tilted
polya-gamma distribution with *b* = 1 and
$c = \sqrt{\mathbb{E}\_{\textbf{z}\_i}(\eta\_{ij}^2)}$, i.e.,
$q(\omega\_{ij}) = PG\\1,\sqrt{\mathbb{E}\_{\textbf{z}\_i}(\eta\_{ij}^2)}\\$.The
terms due to the KL-divergence for *ω*<sub>*i**j*</sub> are given by:

where
𝔼<sub>*ω*<sub>*i**j*</sub></sub>{log (*p*(*ω*<sub>*i**j*</sub>; 1, 0)}
from the prior cancels out with the same term from the entropy of the
variational distribution written as an exponentially tilted polya-gamma
variable with *b* = 1 and *c* = 0. Consequently, the ELBO is now:

so that we see that the terms involving *ω*<sub>*i**j*</sub> cancel out,
and we only need to calculate the normalising constant of the
variational distribution w.r.t. *ω* ∼ *P**G*(1, 0), which has the
solution cosh<sup>−1</sup>(*c*/2).

Finally, with
*η*<sub>*i**j*</sub> = *β*<sub>0*j*</sub> + **z**<sub>*i*</sub><sup>⊤</sup>**γ**<sub>*j*</sub>,
so that $\mathbb{E}\_{\textbf{z}\_i}(\eta\_{ij}) = \tilde{\eta\_{ij}}$,
we have
𝔼<sub>**z**<sub>*i*</sub></sub>(*η*<sub>*i**j*</sub><sup>2</sup>) = *η̃*<sub>*i**j*</sub><sup>2</sup> + tr(**A**<sub>*i*</sub>**γ**<sub>*j*</sub>**γ**<sub>*j*</sub><sup>⊤</sup>).

Plugging in the different results, with the well known results for
$\mathbb{E}\_{\textbf{z}\_i}(\frac{1}{2}\textbf{z}\_i^\top\textbf{z}\_i)$
and 𝔼<sub>**z**<sub>*i*</sub></sub>{log *q*(**z**<sub>*i*</sub>)},
equation becomes:

Finally, we note that for the quadratic model with
*η*<sub>*i**j*</sub> = *β*<sub>0*j*</sub> + **z**<sub>*i*</sub><sup>⊤</sup>**γ**<sub>*j*</sub> − **z**<sub>*i*</sub><sup>⊤</sup>**D**<sub>*j*</sub>**z**<sub>*i*</sub>
the solution to
𝔼<sub>**z**<sub>*i*</sub></sub>(*η*<sub>*i**j*</sub><sup>2</sup>) = var(*η*<sub>*i**j*</sub>) + *η̃*<sub>*i**j*</sub><sup>2</sup>
is the same as for the probit model.
