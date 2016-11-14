---
header-includes:
   - \usepackage{pxfonts}
   - \usepackage{bbm}
---

We model the decision to deworm as
$$ 
Y_{ij}(V^*_{ij}(z), S^*_{ij}(z)) = \mathbbm{1}\{V^*_{ij}(z) + S^*_{ij}(z) > 0\}. 
$$

We further parameterize the distribution of these latent variables to be

$$ 
\begin{pmatrix} V^*_{ij}(control) \\ V^*_{ij}(calendar) \\ V^*_{ij}(bracelet) \\ V^*_{ij}(ink) \end{pmatrix} \sim \mathcal{N} \left( \begin{matrix} \mu_0 \\ \mu_c \\ \mu_b \\ \mu_k \end{matrix},
\begin{matrix} 
  \sigma^2_0 & 0 & 0 & 0 \\
  0          & \sigma^2_c & 0 & 0 \\
  0          & 0 & \sigma^2_b & 0 \\
  0          & 0 & 0 & \sigma^2_k
\end{matrix} \right) 
$$

# Calendar/Bracelet Preference Survey
   
We parameterize the distribution of $D^*_{ij}$, 
$$ 
D^*_{ij} \sim \mathcal{N}(\mu_d, \sigma^2_d), 
$$
such that $\mu_d = \mu_c - \mu_b$ and $\sigma^2 = \sigma^2_c + \sigma^2_b$. Without lost of generality we define $\sigma^2_d = 1$. This allows us to estimate the posterior distribution of $\mu_d$

\begin{equation} f(\mu_d|\mathbf{H}) \label{eqn:choice_posterior} \end{equation}

# Bayesian Analysis

However, we can use the posterior distribution of $D^*_{ij}$ to impute these outcomes.
$$ 
\begin{aligned}
Y_{ij}^{mis,bc} = Y_{ij}(V^*_{ij}(bracelet), S^*_{ij}(calendar)) &= \mathbbm{1}\{V^*_{ij}(bracelet) + S^*_{ij}(calendar) > 0\} \\ 
&= \mathbbm{1}\{V^*_{ij}(calendar) + S^*_{ij}(calendar) + D^*_{ij}> 0\}
\end{aligned}
$$
and
$$ 
\begin{aligned}
Y_{ij}^{mis,cb} = Y_{ij}(V^*_{ij}(calendar), S^*_{ij}(bracelet)) &= \mathbbm{1}\{V^*_{ij}(calendar) + S^*_{ij}(bracelet) > 0\} \\
&= \mathbbm{1}\{V^*_{ij}(bracelet) + S^*_{ij}(bracelet) - D^*_{ij} > 0\}
\end{aligned}
$$

Therefore, the posterior distribution we need in order to estimate \eqref{eqn:social_in_b} and \eqref{eqn:social_in_c} are
$$
f(\mathbf{Y}^{mis}, \mathbf{Y}^{mis,bc}|\mathbf{Y}^{obs}, \mathbf{H}, \mathbf{Z}) = \int_{\theta, \mu_d} f(\mathbf{Y}^{mis}, \mathbf{Y}^{mis,bc}, \theta, \mu_d|\mathbf{Y}^{obs}, \mathbf{H}, \mathbf{Z})\cdot f(\theta, \mu_d| \mathbf{Y}^{obs}, \mathbf{H}, \mathbf{Z})\, \mathrm{d}\theta\mathrm{d}\mu_d
$$
and 
$$
f(\mathbf{Y}^{mis}, \mathbf{Y}^{mis,cb}|\mathbf{Y}^{obs}, \mathbf{H}, \mathbf{Z}) = \int_{\theta, \mu_d} f(\mathbf{Y}^{mis}, \mathbf{Y}^{mis,cb}, \theta, \mu_d|\mathbf{Y}^{obs}, \mathbf{H}, \mathbf{Z})\cdot f(\theta, \mu_d| \mathbf{Y}^{obs}, \mathbf{H}, \mathbf{Z})\, \mathrm{d}\theta\mathrm{d}\mu_d
$$
respectively. Where $\theta = (\mu_{c,total}, \mu_{b,total}, \sigma^2_{c,total}, \sigma^2_{b,total})$ and 
$$
\begin{aligned}
V^*_{ij}(calendar) + S^*_{ij}(calendar) \sim \mathcal{N}(\mu_{c,total}, \sigma^2_{c,total}) \\
V^*_{ij}(bracelet) + S^*_{ij}(bracelet) \sim \mathcal{N}(\mu_{b,total}, \sigma^2_{b,total})
\end{aligned}
$$
