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

### Assignment 1 - Monte Carlo methods  

### Assignmnet 2 - Advanced methods  
