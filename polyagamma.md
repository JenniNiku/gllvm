Let the marginal log-likelihood for a logistic GLLVM be:

\begin{equation}
\begin{aligned}
\mathcal{L}(\Theta) &= \sum \limits ^n_{i=1} \log \displaystyle \int \prod \limits^m_{j=1}p(y_{ij}, \textbf{z}_i)d\textbf{z}_i\\
&= \sum \limits ^n_{i=1} \sum \limits^m_{j=1} \log \mathbb{E}_{\textbf{z}_i}[\exp\{\eta_{ij}\}^{y_{ij}}/\{1+\exp(\eta_{ij})\}h(\textbf{z}_i)],
\end{aligned}
\label{eq1}
\end{equation}

with $\eta_{ij} = \beta_{0j} + \textbf{z}_i^\top\boldsymbol{\gamma}_j$. The above expectation of the first term does not possess a closed-form solution. [Akin to Gibbs sampling of logistic regression](https://arxiv.org/pdf/1205.0310), we apply a polya-gamma data augmentation and additional introduce the auxiliary variable $\omega \sim PG(b,c) = \frac{\exp(-c^2\omega/2)f(\omega;b,0)}{\mathbb{E}\{\exp(-\frac{1}{2}c^2\omega)\}}$, [a polya-gamma random variable](https://stephens999.github.io/fiveMinuteStats/polya_gamma.html), so that the joint log-likelihood instead becomes:

\begin{equation}
\mathcal{L}(\Theta) = \sum \limits ^n_{i=1}\sum \limits^m_{j=1}\log p(y_{ij}\vert\omega_{ij},\textbf{z}_i)+ \log p(\omega;1,0)-\frac{1}{2}\textbf{z}_i^\top\textbf{z}_i,
\label{eq2}
\end{equation}

with corresponding marginal log-likelihood:

\begin{equation}
\mathcal{L}(\Theta) = \sum \limits ^n_{i=1}\sum \limits^m_{j=1}\log \mathbb{E}_{p(\textbf{z}_i)}[\mathbb{E}_{p(\omega_{ij};1,0)}\{p(y_{ij},\omega_{ij},\textbf{z}_i)\}].
\end{equation}

Applying Jensen's inequality, we have the typical VA ELBO:

\begin{equation}
\mathcal{L}_{VA}(\Theta) = \sum \limits ^n_{i=1}\sum \limits^m_{j=1}\mathbb{E}_{\textbf{z}_i,\omega_{ij}}\{\log p(y_{ij}\vert\omega_{ij},\textbf{z}_i)\} + \mathbb{E}_{\omega_{ij}}\{\log p(\omega;1,0)\} -  \mathbb{E}_{\textbf{z}_i}(\frac{1}{2}\textbf{z}_i^\top\textbf{z}_i) - \mathbb{E}_{\omega_{ij}}\{\log q(\omega)\} - \mathbb{E}_{\textbf{z}_i}\{\log q(\textbf{z}_i)\},
\label{eq3}
\end{equation}

where $\mathbb{E}_{\textbf{z}_i}\{\log p(y_{ij}\vert \omega_{ij}, \textbf{z}_i)\} = -\log{2} + (y_{ij} - \frac{1}{2})\mathbb{E}_{\textbf{z}_i}(\eta_{ij}) -\frac{1}{2}\mathbb{E}_{\omega_{ij}}(\omega)\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)$ under the assumption $q(\omega_{ij},\textbf{z}_i) = q(\omega_{ij})q(\textbf{z}_i)$. 

<!-- Note, that we cannot proceed under the assumption $q(\omega_{ij},\textbf{z}_i) = q(\omega_{ij}|\textbf{z}_i)q(\textbf{z}_i)$, as it would leave us with $\mathbb{E}_{\textbf{z}_i}\{\log\cosh(\vert\eta\vert/2)\}$ in the likelihood, which does not possess a closed-form solution. -->

it is straightforward to show that if we formulate $q(\textbf{z}_i\vert \omega_{ij})q(\omega_{ij})$, the first distribution is conditionally normal and thus possesses a closed form solution w.r.t. integration of $\omega_{ij}$. Thus, due to the polya-gamma augmentation, so we can assume (as usual) $q(\textbf{z}_i) = \mathcal{N}(\textbf{a}_i, \textbf{A}_i)$ where $\textbf{a}_i$ and $\textbf{A}_i$ are free variational parameters to estimate. For $q(\omega_{ij})$:

\begin{equation}
\begin{aligned}
\log q(\omega_{ij}) &\propto \mathbb{E}_{\textbf{z}_i}\{\mathcal{L}(\Theta)\}\\
&\propto -\frac{1}{2}\omega\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2) + \log p(\omega;1,0),
\end{aligned}
\end{equation}

which we recognize as the (unnormalized) exponentially tilted polya-gamma distribution with $b = 1$ and $c = \sqrt{\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)}$, i.e., $q(\omega_{ij}) = PG\{1,\sqrt{\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)}\}$.The terms due to the KL-divergence for $\omega_{ij}$ are given by:

\begin{equation}
\begin{aligned}
\mathbb{E}_{\omega_{ij}}\{\log p(\omega_{ij};1,0)\}-\mathbb{E}_{\omega_{ij}}\{\log q(\omega_{ij})\} &= \mathbb{E}_{\omega_{ij}}\{\log p(\omega_{ij};1,0)\}-\mathbb{E}_{\omega_{ij}}[\log\frac{\exp\{-\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)\omega_{ij}/2\}p\{\omega_{ij};1,0\}}{\cosh^{-1}\{\sqrt{\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)}/2\}}]\\
&= \frac{1}{2}\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)\mathbb{E}_{\omega_{ij}}(\omega_{ij})-\log[\cosh\{\sqrt{\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)}/2\}],
\end{aligned}
\end{equation}

where $\mathbb{E}_{\omega_{ij}}\{\log(p(\omega_{ij};1,0)\}$ from the prior cancels out with the same term from the entropy of the variational distribution written as an exponentially tilted polya-gamma variable with $b = 1$ and $c = 0$. Consequently, the ELBO is now:

\begin{equation}
\begin{gathered}
\mathcal{L}_{VA}(\Theta) = -nm\log(2)+\sum \limits ^n_{i=1}\sum \limits^m_{j=1} (y_{ij} - \frac{1}{2}) \mathbb{E}_{\textbf{z}_i}(\eta_{ij}) -\frac{1}{2}\mathbb{E}_{\omega_{ij}}(\omega_{ij})\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2) + \frac{1}{2}\mathbb{E}_{\omega_{ij}}(\omega_{ij})\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2) - \log[\cosh\{\sqrt{\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2)}/2\}] \\+ \sum \limits^n_{i=1}\{\log \det (\textbf{A}_i) - \text{tr}(\textbf{A}_i) - \textbf{a}_i^\top\textbf{a}_i\},
\end{gathered}
\label{eq4}
\end{equation}

so that we see that the terms involving $\omega_{ij}$ cancel out, and we only need to calculate the normalising constant of the variational distribution w.r.t. $\omega \sim PG(1,0)$, which has the solution $\cosh^{-1}(c/2)$.

Finally, with $\eta_{ij} = \beta_{0j} + \textbf{z}_i^\top\boldsymbol{\gamma}_j$, so that $\mathbb{E}_{\textbf{z}_i}(\eta_{ij}) = \tilde{\eta_{ij}}$, we have $\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2) = \tilde{\eta}_{ij}^2+ \text{tr}(\boldsymbol{A}_i\boldsymbol{\gamma}_j\boldsymbol{\gamma}_j^\top)$.

Plugging in the different results, with the well known results for $\mathbb{E}_{\textbf{z}_i}(\frac{1}{2}\textbf{z}_i^\top\textbf{z}_i)$ and $\mathbb{E}_{\textbf{z}_i}\{\log q(\textbf{z}_i)\}$, equation \eqref{eq3} becomes:

\begin{equation}
\begin{aligned}
\mathcal{L}_{VA}(\Theta) =& -nm\log(2)+\sum \limits ^n_{i=1}\sum \limits^m_{j=1} (y_{ij} - \frac{1}{2}) \tilde{\eta_{ij}} - \log\{\cosh(\sqrt{\tilde{\eta_{ij}}^2+\boldsymbol{\gamma}_j^\top\boldsymbol{A}_i\boldsymbol{\gamma}_j}/2)\} + \frac{1}{2}\sum \limits^n_{i=1}\{\log \det (\textbf{A}_i) - \text{tr}(\textbf{A}_i) - \textbf{a}_i^\top\textbf{a}_i\}\\
&= \sum \limits ^n_{i=1}\sum \limits^m_{j=1} (y_{ij} - \frac{1}{2}) \tilde{\eta_{ij}} - \frac{1}{2}\sqrt{\tilde{\eta_{ij}}^2+\boldsymbol{\gamma}_j^\top\boldsymbol{A}_i\boldsymbol{\gamma}_j}-\log\{1+\exp(-\sqrt{\tilde{\eta_{ij}}^2+\boldsymbol{\gamma}_j^\top\boldsymbol{A}_i\boldsymbol{\gamma}_j})\} + \frac{1}{2}\sum \limits^n_{i=1}\{\log \det (\textbf{A}_i) - \text{tr}(\textbf{A}_i) - \textbf{a}_i^\top\textbf{a}_i\}
\end{aligned}
\end{equation}

Finally, we note that for the quadratic model with $\eta_{ij} = \beta_{0j} + \textbf{z}_i^\top\boldsymbol{\gamma}_j - \textbf{z}_i^\top\textbf{D}_j\textbf{z}_i$ the solution to $\mathbb{E}_{\textbf{z}_i}(\eta_{ij}^2) = \text{var}(\eta_{ij})+\tilde{\eta}_{ij}^2$ is the same as for the probit model.