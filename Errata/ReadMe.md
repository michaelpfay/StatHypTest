# List of Corrections to: Statistical Hypothesis Testing in Context
## by Michael P. Fay and Erica H. Brittain


  **p. 206, line 9:**

  The sentence
  "Because $p_{ANOVA}$ is less than $p_{AB}$, $p_{BC}$, and $p_{AC}$, none of the pairwise p-values needs to be adjusted for multiple comparisons."

  Can be more clearly stated as
  "Because $p_{ANOVA}$ is less than $p_{AB}$, $p_{BC}$, and $p_{AC}$, the three unadjusted pairwise p-values equal their adjusted versions." (Thanks to Pekka Pere of Aalto University for pointing this out.)

  **p. 206, line -9:** 
   
   The expression
   "...we test all
   
```math
\begin{pmatrix}
k \\
k-1
\end{pmatrix}
= k-1 
```

null hypotheses that have $k-1$ distributions equal."

should be
"...we test all   

```math
\begin{pmatrix}
k \\
k-1
\end{pmatrix}
= k
```

null hypotheses that have $k-1$ distributions equal." (Thanks to Pekka Pere of Aalto University for catching this error.)

 **p. 207, line 5:** 

 The phrase 
 "Let $\eta(k)$ be the number of tests done...." 
 
 should be 
 "Let $\eta(k)$ be the maximum number of tests done...." (Thanks to Pekka Pere of Aalto University for catching this error.)

   **p. 305, line 6:**

There is a typo in Greenwood's formula, the set in the product operator should exclude $j$ where $d_j=0$. The correct equation is:

$$\widehat{var} \left( \hat{S}(t) \right) =  \left( \hat{S}(t) \right)^2 \prod_{j:t^*_j \leq t \mbox{ and } d_j \neq 0 } \frac{ d_j}{r_j(r_j-d_j)}.$$
