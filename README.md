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

#### 1. Background

##### 1.1 Stock options

Consider the equation for geometric Bornwian motion, as used to model the path of an underlying asset paying proportional dividends at a continuous rate $D_{0}$:
$$dS=(\mu-D_{0})Sdt+\sigma SdW, \qquad (1)$$
where $dW$ is the increment of a Weiner process (drawn from a Normal distribution with mean zero and standard deviation $\sqrt{dt}$); we may then write that
$$dW = \phi\sqrt{dt}, \qquad (2)$$
where $\phi$ is a random variable drawn from a normalised Normal distribution.  
Using (2) and risk neutrality, (1) can be integrated exactly over a timescale $\delta t$ (NOT necessarily small) to yield
$$S(t+\delta t)=S(t)\text{exp}\left[(r-D_{0}-\frac{1}{2}\sigma^{2})\delta t+\sigma\phi\sqrt{\delta t}\right]. \qquad (3)$$
Equation (3) then generates a random path. Since $\delta t$ need not be small, in the case of European options, it is possible to generate a (random) value of S at expiry ($t=T$) in just one step (i.e $\delta t = T$). From this value (say $S(T)$), the payoff can then be easily calculated.  

##### 1.2 European options
Note that because the stock is paying dividends it makes the value of holding a share a little difference since cash dividend payments are made to stock holders. We assume here that all contracts in the portfolio are options, so that \textbf{no} cash payments are received by the owner of the portfolio. Here we may price the options in the portfolio according to the following formula:  
1. Assume that
$$d_{1}=\frac{\text{ln}(S/X)+(r-D_{0}+\sigma^{2}/2)(T-t)}{\sigma\sqrt{T-t}},$$
$$d_{2} = d_{1}-\sigma\sqrt{T-t}.$$
2. A put option P with terminal condition
$$P(S, T) = \text{max}(X-S, 0)$$
has the analytic solution
$$P(S, t) = Xe^{-r(T-t)}N(-d_{2})-Se^{-D_{0}(T-t)}N(-d_{1}).$$
3. A call option C with terminal conditions
$$C(S, T) = \text{max}(S-X,0)$$
has an analytic solution
$$C(S,t)=Se^{-D_{0}(T-t)}N(d_{1})-Xe^{-r(T-t)}N(d_{2}).$$
4. A binary put option BP with terminal conditions
$$BP(S,T) = 1 \qquad \text{if} \qquad S \leq X \qquad \text{or} \qquad 0 \qquad \text{if} \qquad S > X$$
has the analytic solution
$$BP(S,t) = e^{-r(T-t)}N(-d_{2}).$$
5. A binary call option BC with terminal conditions
$$BC(S,T) = 0 \qquad \text{if} \qquad S \leq X \qquad \text{or} \qquad 1 \qquad \text{if} \qquad S > X$$
has the analytic solution
$$BC(S,t) = e^{-r(T-t)}N(d_{2}).$$
6. If there is a payoff at maturity which is equal to the stock price, this is equivalent to the value of a call option with strike ($X=0$) so that
$$C(S, T; X=0) = S$$
and 
$$C(S, t; X=0) = Se^{-D_{0}(T-t)}=e^{-r(T-t)}F_{t,T}$$
which is the discounted futures price $F_{t,T}$.  

If this payoff of a portolio is denoted as $\Pi_{i}(t=T)$ (for the $i^{th}$ simulation), then the value of this payoff at $t=0$ is
$$\Pi_{i}(t=0)=\Pi_{i}(t=T)e^{-rT}. \qquad (4)$$
If N simulations are performed, then we merely average out the $\Pi_{i}(t=0)$ to yield an approximation for the value of the portfolio, i.e
$$\Pi = \frac{\sum_{i=1}{}^{N}\Pi_{i}(t=0)}{N}. \qquad (5)$$

##### 1.2 Path Dependent Options

### Assignmnet 2 - Advanced methods  
