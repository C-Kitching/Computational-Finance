# Computational Finance

A series of computational finance projects based on option pricing. 

## Contents  

### Mini Task 1 - Valuing a simple portfolio  

A trader wants to price the value of a financial contract $\Pi(S, t)$ at time $t=0$ which has the formula  
$$d_{1}=\frac{\text{sinh}\left((S/X)-1\right)+r(T-t)e^{1-(\sigma^{2}/q)}}{\text{exp}[1+\sigma^{2}(T-t)]},$$
$$d_{2}=\frac{\text{sinh}\left((S/X)-1\right)-\sigma\text{sin}(\sigma^{2}-q)\sqrt{T-t}}{\text{exp}[1+\sigma^{2}(T-t)]},$$
$$\Pi(S,t) = Se^{1+\sigma^{2}(T-t)}e^{-r(T-t)}N(d_{1})-\sqrt{X^{1+(r/q)}*S^{1-(r/q)}}e^{-q(T-t)}N(d_{2}),$$
where $T=1$, $X=1500$, $r=0.0319$, $q=0.0207$ and $\sigma = 0.3153$. Here $N(x)$ is the standard normal cummulative distribution function.  
Using the formula and parameters provided, calculate $\Pi$ and output the results to the screen. You should generate 4 columns of data:
1. the value of $S$,
2. the value of $d_{1}$,
3. the value of $d_{2}$,
4. and the value of $\Pi(S, t=0)$.  

Output each of the values when the stock price is
$$S\in [1125,1200,1275,1350,1425,1500,1575,1650,1725,1800,1875].$$

### Mini Task 2 - Valuing a interest rate derivative  

A trader wants to calculate the value of an interest rate derivative contract using a non-standard model. In this particular model, it is assumed that $r(t)$ evolves according to an Ito diffusion process and as such the value of the derivative contract $V(r, t)$ can be found by solving
$$V(r, t) = P(r,t,T)E[g(R_{r,t,T})],$$
where $R_{r,t,T}$ denotes a normal random variable with mean $f(r,t,T)$ and variance $v^{2}(t,T)$. Here $P$ is the value of a pure discount bond paying £1 at time $t=T$ under the risk neutral measure.  
The trader has already solved the SDE under the risk-neutral measure and given you the explicit functions for $f$, $v^{2}(t,T)$ and $P$. They are: 
$$P(r,t,T) = \text{exp}\left[\frac{2}{3}k^{2}(t,T)-\frac{1}{4}n(r,t,T)\right],$$
$$f(r,t,T) = m(r,t,T) - \frac{1}{2}q(t,T),$$
$$v^{2}(t,T) = \frac{\sigma^{2}}{\kappa}(1-e^{-\kappa(T-t)}),$$
where
$$m(r,t,T) = e^{-\kappa(T-t)}r+(1-e^{\kappa(T-t)})\theta,$$
$$n(r,t,T) = r(T-t) - \frac{\theta-r}{2\kappa}(1-e^{-4\kappa(T-t)}),$$
$$k^{2}(t,T) = \frac{\sigma^{2}}{2\kappa^{3}}\left(5e^{-\kappa(T-t)}-3e^{-2\kappa(T-t)}+3\kappa(T-t)-2\right),$$
$$q(t, T) = \frac{\sigma^{2}}{3\kappa^{2}}(1-e^{-\kappa(T-t)})^{5}.$$
The derivative that needs to be priced will be a cash-or-nothing call or put with an interest rate strike price of $X_{r}$. In the case of a call at maturity the holder will recieve £1 is r>X_{r} and nothing otherwise. In the case of a put at maturity the holder will recieve £1 if $r < X_{r}$ and nothing otherwise.  
Then the value of the call option can be written as 
$$V(r,t) = P(r,t,T)(1-N(h)),$$
and a put option is given by
$$V(r,T)=P(r,t,T)N(h),$$
where the value of $h$ must be determined. In order to find h you may use the fact that $R_{r,t,T}$ is normally distributed and $N(x)$ is equal to the probability that $y < x$ when $y$ is drawn from a normal random distribution with mean 0 and variance 1. \it{Hint: you need N(h) to be equivilant to $P(R_{r,t,T} < X).$

### Assignment 1 - Monte Carlo methods  

### Assignmnet 2 - Advanced methods  
