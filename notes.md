
## Reachability of Linear Systems
The computation of the reachability of a dynamic system given the initial conditions is required for the verification of hybrid systems. For this, we will compute tje over-approximation of the reachable sets of uncertain linear systems. An important aspect of this, is the representation of the approximations of the reachable sets. There are multiple methods which use polytopes, oriented hyperrectangles, ellipsoids, orthogonal polyhedra. Here we are going to use orthogonal polyhedra, but using the analysis applied to zonotopes. I'm doing this because, even though the performance is better, implementing the algorithm using zonotopes is harder, but equivalent, so I consider this as a proof-of-concept implementation which could be translated to a zonotope representation with the appropiate methods to transform zonotope-orthogonal polyhedra.

## Approximation of Reachable Sets
Let's start with the following uncertain linear system;

$$ x'(t) = Ax(t) + u(t), ||u(t)|| \leq \mu $$

where A is a $n \times n$ matrix and $||.||$ is the infinity norm on $\mathbb{R}^n$. Given the initial values $I$, if the reachable set at a time t is $\Phi_t(I)$, then the reachable set in the interval $[\underline{t}, \overline{t}]$ can be defined by:

$$ R_{[\underline{t}, \overline{t}]}(I) =  \cup_{t \in [\underline{t}, \overline{t}]} \Phi_t(I)$$

To compute the over-approximation of the reachable set we can decompose $R_{[\underline{t}, \overline{t}]}(I)$ with a time step r.

$$ R_{[0, t]}(I) =  \cup_{i=0}^{i=N-1} R_{[ir, (i+1)r]}(I)$$

So we need to calculate the over approximation set $R_{[ir, (i+1)r]}(I)$. We have

$$ R_{[ir, (i+1)r]}(I) = \Phi_r(R_{[(i-1)r, ir]}(I))$$

## Approximation of $\Phi_r$
Let $y \in \Phi_r$, there is an input u such that

$$ y = e^{rA}x + \int^r_0e^{(r-s)A}u(s)ds$$

so

$$ ||y - e^{rA}x|| \leq \int^r_0e^{(r-s)||A||}\mu ds = \frac{e^{r||A||}-1}{||A||}\mu$$

We can define $\beta = \frac{e^{r||A||}-1}{||A||}\mu$, so $\Phi_r$ is an over approximation of $e^{rA}Z + \square(\beta_r)$

## Approximation of $R_{[0, r]}$

$$R_{[0, r]}(Z) \subseteq (\cup_{t \in [0, r]} e^{rA}Z) + \square(\beta_r)$$

## Mantaining the size of the polyhedra (zonotope)

Because we are doing a Minkowski sum in each iteration the dimension of the Zonotope or the number of vertices of the polyhedra, will keep increasing in each iteration. To avoid this, we are going to use a very simple methods (better methods exist, for example, in [1]). Every time that we have a polyhedraof more than 20 vertices, we are going to substitute it by the smallest parallelogram (the one that covers the polyhedra with the minimum area).

## Reachability Algorithm




```

std::vector< Polygon_2 > Reachability(const Polygon_2& initial_value, double t, 
	double delta_t, const Matrix& a, double mu) {
	int n = t / delta_t;

	double a_norm = a.MaxNorm();
	double alpha_r = (exp(delta_t * a_norm) - 1 - delta_t * a_norm) * Sup(initial_value);
	double beta_r = (exp(delta_t * a_norm) - 1) / a_norm * mu;

	Matrix a_exp = (a * delta_t).Exponential();

	Polygon_2 P0_1 = a_exp * initial_value;

	Matrix a_exp_2 = (a * -delta_t).Exponential();

	Polygon_2 P0_2 = a_exp_2 * initial_value;

	Polygon_2 P0 = MinkowskiSum(Multiply(P0_1, 0.5), Multiply(P0_2, 0.5));

	Polygon_2 Q0 = MinkowskiSum(P0, CreateHyperrectangle(alpha_r + beta_r));

	std::vector< Polygon_2 > res = { Q0 };

	Polygon_2 result = Q0;

	for (int i = 1; i < n; i++) {
		Q0 = ReduceVertices(Q0);

		P0 = a_exp * Q0;
		Q0 = MinkowskiSum(P0, CreateHyperrectangle(beta_r));

		res.push_back(Q0);
	}

	return res;
}
```

## Example

We are going to solve the system:

$$ A = 
\begin{pmatrix}
-1 & -4\\
4 & -1
\end{pmatrix}, \mu = 0.05 $$

with initial value $I = [0.9, 1.1] × [−0.1, 0.1]$.

![alt text](Example_1.png "Title")

## Zonotope

A Zonotope is a set such that

$$ Z = \{ x \in \mathbb{R}^n : x = c + \sum_{i=1}^{i=p}{x_i g_i}, -1 \leq x_i \leq 1  \} $$

where $c, g_i$ are vectors in $\mathbb{R}^n$ called the center and the generators. So Z is defined as $Z=(c, <g_1, g_2, ...>)$.

The Minkowski sum of two zonotopes can be computed as:

$$ Z_1 + Z_2 = (c_1 + c_2, <g_1, g_2, ..., h_1, h_2, ...>)$$

## Reduce Dimension Zonotope 1

The simplest method to compute the vertices of the Zonotope is to calculate all the possible $2^p$ combinations of vertices, and then calculate the convex hull of the vertices. This is computational expensive so it is worth to reduce the dimension of the zonotope.
The method created in [1] is to first sort all the generators by $||Norm_1(g_i) - Norm_{\inf}(g_i)||$ to merge all generators that are almost one dimensional, and then add all the rest of generators. The algorithm is implemented in the method:

```
	Zonotope ReduceVerticesZonotope(const Zonotope& zonotope) {
		if (zonotope.Generators().size() < MAX_VERTICES) return zonotope;

		std::vector<std::vector<double> > generators;

		std::vector<std::vector<double> > old_generators = zonotope.Generators();

		std::sort(std::begin(old_generators),
			std::end(old_generators),
			[](const std::vector<double>& l, const std::vector<double>& r) {
			return L1Norm(l) - LInfNorm(l) < 
				L1Norm(r) - LInfNorm(r);
		});

		int d = old_generators[0].size();
		int p = old_generators.size();
		int r;
		if (d == 2) {
			r = 4;
		}
		else {
			r = 2;
		}
		int m = p - d * (r - 1);
		for (int i = 0; i < d; i++) {
			std::vector<double> v(d, 0.0);
			for (int j = 0; j < m; j++) {
				v[i] += abs(old_generators[j][i]);
			}
			generators.push_back(v);
		}

		for (int i = m; i < p; i++) {
			generators.push_back(old_generators[i]);
		}
		
		return Zonotope(zonotope.Center(), generators);
	}
```

## Equations of the Car

$$ \dot{x} = \dot{x_1} = s \cdot cos(x_3 + x_5)$$
$$ \dot{y} = \dot{x_2} = s \cdot sin(x_3 + x_5)$$

$$ \dot{\theta} = x_4$$

$$ \ddot{\theta} = \dot{x_4} = -\frac{c_1}{s}x_4 - c_2x_5 + c_3u$$

$$ \dot{\phi} = \dot{x_5} = (-1 - \frac{c_4}{s^2})x_4 - \frac{c_5}{s}x_5 + \frac{c_6}{s}u$$

## Linearization of the Equations

The previous system of equations is non-linear, but we need to use a linear system to solve it and so we can separate the computation of the reachable sets of the state and input dependent part.

The dynamics of the autonomous car is given by:

$$ \dot{x} = f(x, u)$$

$$ ||u_{\inf}|| < \delta $$

We can approximate this by a linear system by keeping the first order of the Taylor expansion:

$$ \dot{x} \approx f(x_i^*,u^*) + \frac{\partial f(x, u)}{\partial x}\bigg\rvert_{x=x_i^*,u=u^*} \Delta x + \frac{\partial f(x, u)}{\partial u}\bigg\rvert_{x=x_i^*,u=u^*} \Delta u  + ...$$

where we can define:

$$ A_i = \frac{\partial f(x, u)}{\partial x}\bigg\rvert_{x=x_i^*,u=u^*}$$

$$B_i = \frac{\partial f(x, u)}{\partial u}\bigg\rvert_{x=x_i^*,u=u^*}$$

and $x_i^*=mid(X_i), u^*=0$

## Reachable Set Computations

There are two types of reachable sets, the state-dependent and input-dependent sets.

| Reachable set        | x(0)           | $f(x^*, u^*)$  | u |
| ------------- |:-------------:| -----:| -----:|
| R      | $\neq 0$ | $\neq 0$ | $\neq 0$ |
| $\bar{R}$      | $= 0$      |   $= 0$ | $\neq 0$ |
| $\hat{R}$ | $\neq 0$      |    $= 0$ | $= 0$ |
| $\check{R}$ | $= 0$      |    $\neq 0$ | $= 0$ |

### State dependent Reachable set
$\hat{R}$ is computed with the method for linear systems described above.

$\check{R}$ is computed with:

$$ \check{R}_i(T) = \int_0^T{a^{A_i(t-\tau)}d\tau f(x_i^*, u*)} = A_i^{-1}(e^{A_i T} - I)f(x_i^*, u*)$$

Ans $\check{R}([0, T])$ is approcimated by:

$$ \check{R}([0, T]) = \frac{t}{T} \check{R}_i(T) + \square (\eta_i)$$

$$ \eta_i = || A_i^{-1}[e^{A_i T} - I - A_iT - 0.375(A_iT)^2]f(x_i^*, u*)||_\infin$$

### Input dependent Reachable set
$\bar{R}$ is computed by using u as a disturbance.

$$ \bar{R}(kT) = A_i^{-1}(e^{A_iT} - I)B_iu(k) + A_i^{-2}(e^{A_iT} - I - A_iT)B_i\dot{u}(k)$$

## Separation of Linear System

If we use the equation as we have defined them earlier, we will notice the the matrix A is non-invertible, but we are making use of the inverse in the calculation of the reachable sets. 

To solve this, and considering the equations are coupled (($x_4, x_5$) -> ($x_3$)-> ($x_1, x_2$)), we can solve the systems separately and create three Marokov Chains.


# Markov Chain

We have calculated the reachable sets $R_i(T)$ and $R_i([0, T])$, and now we need to calculate the transition probabilities of discrete time Markov Chains that are use to model the continuous linear dynamics of the model.

We are going to define the probability that the state is in a cell $X_i$ at a point of time or time interval. The transitions probabilities of the state dependant Markov chains are:

$$ \Phi_{ij}^*(T) = \frac{V((\hat{R}_j(T) \oplus \check{R}_j(T)) \cap X_i)}{V(\hat{R}_j(T) \oplus\check{R}_j(T))} $$

$$ \Phi_{ij}([0, T]) = \frac{V(R_j([0, T])  \cap X_i)}{V(R_j([0, T]))} $$

where V() is the volume of the polytope.

The transition Matrix for the input-dependant solution is:

$$ \bar{\Phi}_{ij}(T) = \frac{V((R_j(T)  \oplus X_j) \cap X_i)}{V(X_j)} $$

Then with the transitions matrixes we can compute the probabilities.

$$ p_i((k+1)T) = \Phi_{il}^*(T) \bar{\Phi}_{lj}(T) p_j(kT)$$

$$ p_i([kT, (k+1)T]) = \Phi_{ij}([0, T])p_j(kT)$$

## Separation of the systems of differential equations

The method explained here can't be applied directly to the system described here for two main reasons:

- The matrix A is non-invertible, so we can't compute the reachability sets directly.
- There are 5 variables, so the discretization with n elements would give an algorithm and memory complexity O(n^5) which is computationally very expensive.

The good thing is that this system of differential equations can be separated in 3 independent systems as $(x_4, x_5)$ -> $(x_3)$ -> $(x_1, x_2)$. 

We only need to calculate the reachable set for $(x_4, x_5)$. Once we have calculated this, we can compute the Markov chains probability matrices and the probabilities for each cell.

If we have the probability for each state $(x_4, x_5)$, we can calculate the probability of the other state cells $(x_3)$, $(x_1, x_2)$, by propagating the differential equations, and because they are independent, then we can merge the probabilities of the states as:

$$p(x_1, x_2, x_3, x_4, x_5) = p(x_1, x_2)\cdotp(x_3)\cdot p(x_4, x_5)$$


## Results


