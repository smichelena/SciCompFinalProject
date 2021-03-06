### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ bfef8ade-0b16-4a28-86d7-7e91876519c9
begin
	using ShortCodes, ExtendableGrids, VoronoiFVM, PlutoVista
	using Plots, GridVisualize, PlutoUI, LaTeXStrings, PyPlot
	using PyCall, CSV, DataFrames
	GridVisualize.default_plotter!(Plots);
	pyplot();
	#TableOfContents();
end

# ╔═╡ c3c70bce-7f23-477b-8f6d-d51b42ed6b2e
md"""
# Bidomain problem
"""

# ╔═╡ 6c8a63a7-2f5f-4e10-aca5-06d58d0eacd1
md"""
# Table of Contents
1. Introduction
1. Bidomain Problem
1. Biodomain Problem Modeled as a Partial Differential Equation
1. Discretization Stategy
1. Results
1. References


# Introduction
This document shall document and implement our strategy for solving the "bidomain problem". It introduces the problem at hand, illustrates our discretization strategy, displays our implementation of this strategy, and visualizes the resulting solutions to the problem.

The solution to the problem, once set up, is calculated by use of the VoronoiFVM.jl Julia package [1]. This package is available at **https://github.com/j-fu/VoronoiFVM.jl**.

The bidomain problem was chosen as the project examination topic by our group as part of the Scientific Computing course at the Teschnische Universität Berlin. This work was performed during the Winter Semester 2021/2022.

# Bidomain problem
In summary, the bidomain problem is a system of partial differential equations modeling the propagation of electric signals throughout cardiac tissue[2]. Before going into more detail, it is convenient to provide some context on the lower fidelity, but simpler, "monodomain problem". 

The monodomain problem models cardiac tissue's electircal properties by treating it as an intracellular region separated from a extracellular region (the electrical ground) by a membrane. This approach neglects current flow external to the cells [2].

The bidomain model then rectifies this issue by treating the intracellular and extracellular regions with their own current flows, linked by transmembrane potential. This allows the regions to be modeled with separate boundary conditions and for the application of external stimuli [2].

This document will use the formulation of the biodomain problem provided in Ether and Bourgalt, *Semi-Implicit Time-Discretization Schemes For The Bidomain Model* [3].
"""

# ╔═╡ 9653fb91-a348-4c4e-9892-020211969393
md"""
# Biodomain Problem Modeled as a Partial Differential Equation
The Bidomain problem for membrane-based models is given by Ethier & Bourgalt [3] as
```math
\newcommand{\math}[1]{\begin{align}#1\end{align}}
\newcommand{\pth}[1]{\left(#1\right)}
\newcommand{\brc}[1]{\left\{#1\right\}}
\newcommand{\bra}[1]{\left[#1\right]}
\newcommand{\dd}{\text{d}}
\newcommand{\eps}{\varepsilon}
\math{
	\frac{\partial u}{\partial t} &= \frac{1}{\eps}f(u,v)+
	\nabla\cdot(\sigma_i\nabla u)+\nabla\cdot(\sigma_i\nabla u_e)\\
	0&=\nabla\cdot(\sigma_i\nabla u+(\sigma_i +\sigma_e)\nabla u_e)\\
	\frac{\partial v}{\partial t} &= \eps g(u,v).
}
```
where [3]:
* ``u = u_i - u_e`` is the transmembrane potential
* ``\sigma_i`` and ``\sigma_e`` are second order tensors representing the intracellular and extracellular tissue's electrical conductivity in each spatial direction
* ``v`` is a lumped ionic variable
* ``\eps`` is a parameter linked to the ratio between the repolarisation rate and the tissue excitation rate
"""

# ╔═╡ 41ad9bd9-6cea-4b44-a3af-297b21db9f67
md"""
Which can also be written in vector form as
```math
\newcommand{\u}{\vec{u}}
\newcommand{\r}{\vec{r}}
\newcommand{\s}{\vec{s}}
\newcommand{\j}{\vec{j}}
\newcommand{\f}{\vec{f}}
\math{
	\partial_t \s(\u) + \nabla \cdot \j(\u) + \r(\u) = \f
}
```
where
```math
\newcommand{\mtr}[1]{\begin{bmatrix}#1\end{bmatrix}}
\math{
	\u &= \mtr{u\\u_e\\v}
	\\\s(\u) &= \mtr{u\\0\\v}\\
\j(\u) &= \mtr{-\sigma_i(\nabla u+\nabla u_e)\\
\sigma_i\nabla u + (\sigma_i+\sigma_e)\nabla u_e
\\0}\\
\r(\u) &= \mtr{-f(u,v)/\eps\\0\\ -\eps g(u,v)}\\
\f &= \vec{0}.
}
```
"""

# ╔═╡ 5c0fa52c-aad8-46c3-a660-4ac85b3904ea
md"""
With $f$ and $g$ defined as
```math
\math{
f(u,v) &= u - \frac{u^3}{3} - v\\
g(u,v) &= u + \beta -\varphi v.
}
```
"""

# ╔═╡ b08f1ba2-ff89-4d25-84cf-ef100cde2d91
md"""
The initial condition for this problem in 1D will be set so that
```math
\math{
	u=u_{0,1},\;v&=v_0 \Rightarrow f(u_{0,1},v_0) = g(u_{0,1},v_0) = 0 \;\;\forall x\notin [0,L/20]
	\\
	u=u_{0,2} = 2,\;v&=v_0 \;\;\forall x\in [0,L/20]\\
	u_e = u_{e,0} &=0\;\;\forall x.
}
```
"""

# ╔═╡ eb03ae91-d8bf-43cd-87b4-9a925b5051f4
md"""
We can get the exact solution to $f=g=0$:
```math
\math{
	0 &= f(u_{0,1},v_0) =  u_{0,1} - \frac{u_{0,1}^3}{3} - v_0 \Rightarrow v_0 = u_{0,1} - \frac{u_{0,1}^3}{3}\\
	0&=  g(u_{0,1},v_0) = u_{0,1} + \beta - \varphi v_0 \Rightarrow  0 = u_{0,1} + \beta -\varphi\pth{u_{0,1} - \frac{u_{0,1}^3}{3}}\Rightarrow\\
	0 &=  u_{0,1}^3\frac{\varphi}{3} +u_{0,1}\pth{1-\varphi} +\beta
}
```
which has three solutions since it is a polynomial of degree 3. Two of them are complex, whilst the third one is real. We find the solution by using Newton's method since the expression for the solution of a 3rd degree polynomial is very verbose.
"""

# ╔═╡ 60a981c1-5e6b-459b-b7b1-dab2b94e1881
md"""
### Finite Volume Discretization
When doing the finite volume discretization we need to split the polygonal domain $\Omega$ into finite volumes $\omega_k$ such that
```math
\math{
	\bar\Omega &= \bigcup_{k\in N} \bar\omega_k\\
	N &= \brc{1,\dots,n_c}\\
	\partial\Omega &= \bigcup_{m\in G}\Gamma_m\\
	G &= \brc{1,\dots,n_e}\\
	\s_{kl} &= \bar\omega_k\cup\bar\omega_l
}
```
where $n_c$ is the number of control volumes, $n_e$ the number of edges of the polygonal domain. We have that $|s_{kl}|>0$ makes $\omega_k$,$\omega_l$ neighbours. We define the part of $\omega_k$ that is on the boundary of $\Omega$ as
```math
\math{
	\gamma_{km} = \partial\omega_k \cup \Gamma_m
}
```
"""

# ╔═╡ 3ee2e07c-3463-470b-a478-bc7e2464389d
md"""
#### Discretization of equation 1
Integrate over a control volume $\omega_k$:
```math
\begin{equation}
   \int_{\omega_k}\frac{\partial u}{\partial t}=\int_{\omega_k}\frac{1}{\eps}f(u,v)d\omega +
       \int_{\omega_k}\nabla \cdot (\sigma_i \nabla u)d\omega + \int_{\omega_k}\nabla \cdot (\sigma_i \nabla u_e)d\omega
\end{equation}
```
Apply Gauss' thereom to the two divergence terms:
```math
\begin{equation}
    \int_{\omega_k}\frac{\partial u}{\partial t} =\int_{\omega_k}\frac{1}{\eps}f(u,v)d\omega +
    \int_{\partial\omega_k}\sigma_i \nabla u \cdot \vec{n}ds +
    \int_{\partial\omega_k}\sigma_i \nabla u_e \cdot \vec{n}ds
\end{equation}
```
Convert the surface integrals over the edge of the control volume into a sum of an integral over each side, as our Voronoi cell control volumes are polygons. Additionally, introduce a separate term for the surface integral over the specific case of sides that are boundary conditions.
```math
\begin{align*}
    \int_{\omega_k}\frac{\partial u}{\partial t} &=\int_{\omega_k}\frac{1}{\eps}f(u,v)d\omega \\ 
    &+ \sum_{l \in N_k}\int_{s_{kl}} \sigma_i \nabla u \cdot \vec{n}_{kl}ds +
    \sum_{m \in B_k}\int_{\gamma_{kl}} \sigma_i \nabla u \cdot \vec{n}_{m}ds \\
    &+ \sum_{l \in N_k}\int_{s_{kl}} \sigma_i \nabla u_e \cdot \vec{n}_{kl}ds +
    \sum_{m \in B_k}\int_{\gamma_{kl}} \sigma_i \nabla u_e \cdot \vec{n}_{m}ds
\end{align*}
```
Exploiting the admissibility condition, approximate the dot products using finite differences:
```math
\begin{equation}
    \sigma_i \nabla u \cdot \vec{n} = \sigma_i \frac{u_k - u_l}{|x_k - x_l|}
    \sigma_i \nabla u_e \cdot \vec{n} = \sigma_i \frac{u_{e_k} - u_{e_l}}{|x_k - x_l|}
\end{equation}
```
Substitute this approximation in and replace the integrals with a simple multiplication by the length of the cell border. Additionally expanding $f(u,v)$:
```math
\begin{align*}
|\omega_k|\frac{\partial u_k}{\partial t} &= \int_{w_k}\frac{1}{\eps}(u_k - \frac{u_k^3}{3} - v_k)d\omega\\
&+ \sum_{l \in N_k} |s_{kl}| \sigma_i  \frac{u_k - u_l}{|x_k - x_l|}  + 
\sum_{m \in B_k} |\gamma_{km}| \sigma_i  \frac{u_k - u_l}{|x_k - x_l|}   \\
&+ \sum_{l \in N_k} |s_{kl}| \sigma_i  \frac{u_{e_k} - u_{e_l}}{|x_k - x_l|}  + 
\sum_{m \in B_k} |\gamma_{km}| \sigma_i  \frac{u_{e_k} - u_{e_l}}{|x_k - x_l|}  
\end{align*}
```
Combine terms, and note the first integral simply becomes a multiplication by $|\omega_k|$:
```math
\begin{align*}
|\omega_k|\frac{\partial u_k}{\partial t} &= \frac{|\omega_k|}{\eps}\pth{u - \frac{u^3}{3} - v}\\
&+ \sum_{l \in N_k} |s_{kl}| \sigma_i  \frac{u_k - u_l + u_{e_k} - u_{e_l}}{|x_k - x_l|} \\
&+ \sum_{m \in B_k} |\gamma_{km}| \sigma_i  \frac{u_k - u_l + u_{e_k} - u_{e_l}}{|x_k - x_l|}  
\end{align*}
```
Move the boundary terms to the other side and organize so it is on the right:
```math
\begin{align}
|\omega_k|\frac{\partial u_k}{\partial t} + \sum_{l \in N_k} |s_{kl}| \sigma_i \frac{u_k - u_l + u_{e_k} - u_{e_l}}{|x_k - x_l|} 
&+ \frac{|\omega_k|}{\eps}\pth{u_k - \frac{u_k^3}{3} - v_k}
=\\
-\sum_{m \in B_k} |\gamma_{km}| \sigma_i  \frac{u_k - u_l + u_{e_k} - u_{e_l}}{|x_k - x_l|} 
\end{align}
```
Then discretizing in time yields:
```math
\begin{align}
\frac{|\omega_k|}{\Delta t}(u_k^n - u_k^{n-1}) + \sum_{l \in N_k} |s_{kl}| \sigma_i \frac{u_k^{\theta} - u_l^{\theta} + u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|} 
&+ \frac{|\omega_k|}{\eps}\pth{u_k^{\theta} - \frac{(u_k^{\theta})^3}{3} - v_k^{\theta}}
=\\
-\sum_{m \in B_k} |\gamma_{km}| \sigma_i  \frac{u_k^{\theta} - u_l^{\theta} + u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|} 
\end{align}
```
Where $u_k^{\theta} = \theta u_k^n + (1+\theta) u_k^{n-1}$
"""

# ╔═╡ 2a0bc645-8dec-454e-9eb8-623c7e95a20c
md"""
#### Discretization of equation 2
Distrbute the divergence and integrate over a volume $\omega_k$:

```math
\begin{align*}
0 = \int_{\omega_k} \nabla \cdot \pth{\sigma_i \nabla u}d\omega + \int_{\omega_k}\nabla \cdot \pth{\sigma_i + \sigma_e}\nabla u_e d\omega
\end{align*}
```
Apply Gauss' theorem:
```math
\begin{align*}
0 = \int_{\partial\omega_k} \pth{\sigma_i \nabla u}\cdot \vec{n}ds + \int_{\partial \omega_k}\pth{\sigma_i + \sigma_e}\nabla u_e \cdot \vec{n} ds
\end{align*}
```
Again convert these terms to a sum of the integral over each side of the volume, and add terms for the boundary conditions:
```math
\begin{align*}
    0 &= \sum_{l \in N_k}\int_{s_{kl}} \sigma_i \nabla u \cdot \vec{n}_{kl}ds +
    \sum_{m \in B_k}\int_{\gamma_{kl}} \sigma_i \nabla u \cdot \vec{n}_{m}ds \\
    &+ \sum_{l \in N_k}\int_{s_{kl}} \pth{\sigma_i + \sigma_e} \nabla u_e \cdot \vec{n}_{kl}ds +
    \sum_{m \in B_k}\int_{\gamma_{kl}} \pth{\sigma_i + \sigma_e}\nabla u_e \cdot \vec{n}_{m}ds
\end{align*}
```
Again approximate the dot products as finite differences, and replace the integrals by the length of the cell border:
```math
\begin{align*}
0 &= \sum_{l \in N_k} |s_{kl}| \sigma_i \pth{ \frac{u_k - u_l}{|x_k - x_l|} } + 
\sum_{m \in B_k} |\gamma_{km}| \sigma_i \pth{ \frac{u_k - u_l}{|x_k - x_l|} }  \\
&+ \sum_{l \in N_k} |s_{kl}| \pth{\sigma_i+\sigma_e} \pth{ \frac{u_{e_k} - u_{e_l}}{|x_k - x_l|} } + 
\sum_{m \in B_k} |\gamma_{km}| \pth{\sigma_i+\sigma_e} \pth{ \frac{u_{e_k} - u_{e_l}}{|x_k - x_l|} }
\end{align*}
```
Move the boundary conditions to the opposite side:
```math
\begin{align*}
\sum_{l \in N_k} |s_{kl}| \sigma_i \frac{u_k - u_l}{|x_k - x_l|}
+ \sum_{l \in N_k} |s_{kl}| \pth{\sigma_i+\sigma_e}  \frac{u_{e_k} - u_{e_l}}{|x_k - x_l|}  = \\
-\sum_{m \in B_k} |\gamma_{km}| \sigma_i \frac{u_k - u_l}{|x_k - x_l|}  
-\sum_{m \in B_k} |\gamma_{km}| (\sigma_i+\sigma_e)  \frac{u_{e_k} - u_{e_l}}{|x_k - x_l|}
\end{align*}
```
Then discretizing in time yields:
```math
\begin{align}
\sum_{l \in N_k} |s_{kl}| \sigma_i \frac{u_k^{\theta} - u_l^{\theta}}{|x_k - x_l|}
+ \sum_{l \in N_k} |s_{kl}| \pth{\sigma_i+\sigma_e}  \frac{u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|}  = \\
-\sum_{m \in B_k} |\gamma_{km}| \sigma_i \frac{u_k^{\theta} - u_l^{\theta}}{|x_k - x_l|}  
-\sum_{m \in B_k} |\gamma_{km}| (\sigma_i+\sigma_e)  \frac{u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|}
\end{align}
```
"""

# ╔═╡ 2cb6b44b-6236-40e3-a945-4bb8df3cae1e
md"""
#### Discretization of equation 3
With 
```math
\begin{align*}
\frac{\partial v}{\partial t} - \eps(u + \beta - \varphi v) = 0
\end{align*}
```
We take the integral as before with respect to the volume $\omega_k$:
```math
\begin{align*}
\int_{\omega_k}\frac{\partial v}{\partial t} d\omega - \int_{\omega_k} \eps(u + \beta - \varphi v)d\omega = 0
\end{align*}
```
The integral is simply a multiplication by the area of the volume:
```math
\begin{align*}
|\omega_k|\frac{\partial v_k}{\partial t} - \eps |\omega_k| (u_k + \beta - \varphi v_k) = 0
\end{align*}
```
Then discretizing in time yields:
```math
\begin{align*}
\frac{|\omega_k|}{\Delta t}(v_k^n - v_k^{n-1}) - \eps |\omega_k| (u_k^{\theta} + \beta - \varphi v_k^{\theta}) = 0
\end{align*}
```

This concludes the space discretization.
"""

# ╔═╡ d7aeab82-eab7-44f1-8586-4504287a4dfc
md"""
### Time discretization

There are a few possibilities for the time discretization:

"""

# ╔═╡ 042f94d9-04d0-4a03-b3b6-375f8b22953b
md"""
#### Implicit Euler
With an implicit euler we approximate $\frac{\partial v}{\partial t}$ by a finite difference between the current time step and the last time step. That is, we set $\theta = 0$

This provides an easy implementation as all we need to implement this approximation is a storage of the previous $\mathbf{u}$. We can then simply solve the resulting system to find the current $\mathbf{u}$.
"""

# ╔═╡ ac978664-28a5-40bb-8a99-4efbf3c6c6fe
md"""
#### Explicit Euler
The explicit or forward Euler yields difficulties when implementing because the time step size must then be much smaller than the grid resultion. In fact, According to the paper by Ethier and Bourgalt, we must have $\Delta t \in O(\min(\eps/L_f, \frac{m_i^3}{M_i^4}h^2))$ where $L_f$ is a Lipschitz constant for $f$ and where $m_i$ and $m_e$ are the 1-ellipticity constants for $u$ and $u_e$.

Note: This analysis is done with a system discretized using FEM instead of FVM, however, its reasonable to assume the stability results are also valid here, as the system properties remain the same irregardless the discretization method applied.
"""

# ╔═╡ 63d685c1-dcd1-454f-8e56-f697e7aca8fb
md"""
### Discretization summary
In total, for the discretization of the problem, we are left with a nonlinear system of equations that must be solved at each time step:

```math
\begin{align}
\frac{|\omega_k|}{\Delta t}(u_k^n - u_k^{n-1}) + \sum_{l \in N_k} |s_{kl}| \sigma_i \frac{u_k^{\theta} - u_l^{\theta} + u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|} 
&+ \frac{|\omega_k|}{\eps}\pth{u_k^{\theta} - \frac{(u_k^{\theta})^3}{3} - v_k^{\theta}}
=\\
-\sum_{m \in B_k} |\gamma_{km}| \sigma_i  \frac{u_k^{\theta} - u_l^{\theta} + u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|} 
\end{align}
```
```math
\begin{align}
\sum_{l \in N_k} |s_{kl}| \sigma_i \frac{u_k^{\theta} - u_l^{\theta}}{|x_k - x_l|}
+ \sum_{l \in N_k} |s_{kl}| \pth{\sigma_i+\sigma_e}  \frac{u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|}  = \\
-\sum_{m \in B_k} |\gamma_{km}| \sigma_i \frac{u_k^{\theta} - u_l^{\theta}}{|x_k - x_l|}  
-\sum_{m \in B_k} |\gamma_{km}| (\sigma_i+\sigma_e)  \frac{u_{e_k}^{\theta} - u_{e_l}^{\theta}}{|x_k - x_l|}
\end{align}
```
```math
\begin{align}
\frac{|\omega_k|}{\Delta t}(v_k^n - v_k^{n-1}) - \eps |\omega_k| (u_k^{\theta} + \beta - \varphi v_k^{\theta}) = 0
\end{align}
```
This yields a system of $\sum_{k \in N}\sum_{l \in N_k}6$ nonlinear equations to be solved in each time step.

"""

# ╔═╡ b79e3858-2208-4cf9-9834-c9bc1c8e3486
md"""
# Possible solution methods

Possible solution methods for our aquired nonlinear system could be fixpoint iteration, which gives a large area of convergence, but it could be really slow. Other methods include Newton iteration which would need a better initial guess for convergence, but the convergence near the solution is quadratic. 

To simplify Newton iteration you can use dual numbers for automatic differentiation, as the VoronoiFVM package does.

Additionally, for performance improvements, you can use a damped Newton scheme or even an embedded newton scheme if the problem is parameter dependant.

"""

# ╔═╡ 3e1c814e-f69d-4839-8ca1-ed03952188a8
md"""
# Discrete solution using VoronoiFVM
### 1D Problem
"""

# ╔═╡ 71943829-1f12-46ce-8347-99a4aebac3d4
begin
	β = 1.0
	γ = 0.5
	ε = 0.1
	σᵢ_normal = 1.0
	σₑ_normal = 1.0
	σᵢ_anisotropic = 25*[0.263 0; 0 0.0263]
	σₑ_anisotropic = 25*[0.263 0; 0 0.1087]
end;

# ╔═╡ 63372ffa-406b-4eef-a9b2-394326cd35a6
begin
	L = 70
	T = 30
end;

# ╔═╡ d3b4bdc4-f17f-4a83-9ddc-bfdd1ea7f26c
f(u,v) = u - u^3/3 - v

# ╔═╡ d37259c3-f3f7-443c-87fb-497f700e6bac
g(u,v) = u + β - γ*v

# ╔═╡ 78326e2d-4e56-4b4a-8b4f-6f31662b3f8e
"""
	fg_zeros()

Returns the approximate solution to f=g=0.
"""
function fg_zeros()
	p = 3(1-γ)/γ; q = 3β/γ;
	h(u) = u^3 + p*u + q; ḣ(u) = 3u^2 + p
	uₙ = 1; hₙ = 1
	cnt = 0
	while abs(hₙ) > 1e-15
		uₙ = uₙ - h(uₙ)/ḣ(uₙ)
		hₙ = h(uₙ)
	end
	vₙ = uₙ - uₙ^3/3
	uₙ, vₙ
end

# ╔═╡ b9adb394-034e-4ec5-a539-2a0e850fe1b5
u₀, v₀ = fg_zeros()

# ╔═╡ 7a7229fc-771c-4f87-a705-a29ddc546ac1
"""
	create_grid(N, dim=1)

Creates a simplex grid in dimension `dim` with `N` grid points in one
dimension with length `L`.
"""
function create_grid(N, dim=1)
	if dim == 2
		X = LinRange(0,L,N[1])
		Y = LinRange(0,L,N[2])
		simplexgrid(X,Y)
	else
		X = LinRange(0,L,N)
		simplexgrid(X)
	end
end

# ╔═╡ 761e9442-62ab-41da-86e3-cc5ef768dd49
md"""
	Functions for getting initial conditions
"""

# ╔═╡ 47c78a62-d0a1-4120-9303-35a9d7d20ee9
begin
	function u⃗₀(x)
		u = u₀; v = v₀
		if 0 ≤ x ≤ L/20
			u = 2
		end
		uₑ = 0
		return [u, uₑ, v]
	end
	function u⃗₀(x,y)
		u⃗₀(x)
	end
	function u⃗₀_2D(x,y)
		u = u₀; v = v₀
		if 0 ≤ x ≤ 3.5 && 0 ≤ y ≤ 70
			u = 2
		end
		if 31 ≤ x ≤ 39 && 0 ≤ y ≤ 35
			v = 2
		end
		uₑ = 0
		return [u, uₑ, v]
	end
end

# ╔═╡ dea5cc0a-a112-494a-85ba-0fcf15d6c2e6
let
	X=0:0.1:L
	initial = u⃗₀.(X)
	init = [lst[k] for lst in initial, k=1:3]'
	u₀ = init[1,:]; uₑ₀ = init[2,:]; v₀ = init[3,:];
	lw = 3
	p = Plots.plot(title="Initial conditions",
		legend=:topright, legendfontsize=10,
		xlabel=L"x", ylabel=L"y",
		xtickfontsize=10, ytickfontsize=10, 
		xguidefontsize=15, yguidefontsize=15,
		yguidefontrotation=-90);
	Plots.plot!(p,X,u₀,label=L"u_0",lw=lw, color="Blue")
	Plots.plot!(p,X,uₑ₀,label=L"u_{e,0}",lw=lw, color="Red")
	Plots.plot!(p,X,v₀,label=L"v_0",lw=lw, color="Green")
end

# ╔═╡ f4f6357b-3d1b-4730-909c-8a62ae8c284e
md"""
	function for getting the problem physics depending on the setup
"""

# ╔═╡ 56bec1ef-6b16-4a37-9566-f5fb4e6a2190
md"""
	Function for getting the anisotropic flux
"""

# ╔═╡ 1fc8a5f0-fa99-44c6-a19a-3334e20eec15
function get_σₛ(σᵢ,σₑ,edge)
	n1 = edge.coord[:,edge.node[1]]
	n2 = edge.coord[:,edge.node[2]]
	x = abs(n1[1] - n2[1]); y = abs(n1[2] - n2[2])
	θ = atan(y/x)
	l = [1 1]; r = [cos(θ);sin(θ)]
	σᵢ_θ = (l*σᵢ*r)[1]; σₑ_θ = (l*σₑ*r)[1]
	σᵢ_θ, σₑ_θ
end

# ╔═╡ 6cb8de00-62a6-46b0-a367-21a282e4ff89
function create_physics(σᵢ_orig, σₑ_orig; initial2D=false, anisotropic=false)
	physics = VoronoiFVM.Physics(
		storage = function(y,u,node)
			y[1] = u[1]
			y[2] = 0
			y[3] = u[3]
		end,
		flux = function(y,u,edge)
			σᵢ, σₑ = anisotropic ? get_σₛ(σᵢ_orig,σₑ_orig,edge) : (σᵢ_orig,σₑ_orig)
			y[1] = -σᵢ*(u[1,2]-u[1,1] + u[2,2]-u[2,1])
			y[2] = σᵢ*(u[1,2]-u[1,1]) + (σᵢ+σₑ)*(u[2,2]-u[2,1])
			y[3] = 0
		end,
		reaction = function(y,u,node)
			y[1] = -f(u[1],u[3])/ε
			y[2] = 0
			y[3] = -ε*g(u[1],u[3])
		end,
		breaction = function(y,u,node)
			if node.coord[:,node.index] == [0, 0]
				y[2] = 0
			end
		end,
	)
end

# ╔═╡ d5315aef-978d-447f-85dc-a3becdc33078
function save_initial(init, filename)
	mkpath("../csv")
	df = DataFrame(init, :auto)
	CSV.write("../csv/$filename", df)
end

# ╔═╡ f19b980c-7d6c-41c6-99af-8475f7aa72db
function get_initial_data(filename)
	Array(DataFrame(CSV.File("../csv/$filename")))
end

# ╔═╡ 6ef697e9-a4b8-446e-bbb8-b577cea0161d
function initial_cond(sys, xgrid, Δt, Tinit_solve, 
		dim, initial2D, use_csv, anisotropic)
	init = unknowns(sys)
	U = unknowns(sys)
	filename="initial_2D"*(anisotropic ? "_anisotropic" : "")*".csv"
	if dim==2 && initial2D
		if use_csv
			init = get_initial_data(filename)
		else
			inival = map(u⃗₀_2D, xgrid)
			init .= [tuple[k] for tuple in inival, k in 1:3]'
			for t ∈ 0:Δt:Tinit_solve
				solve!(U, init, sys; tstep=Δt)
				init .= U
			end
			save_initial(init, filename)
		end
	else
		inival = map(u⃗₀, xgrid)
		init .= [tuple[k] for tuple in inival, k in 1:3]'
	end
	init
end

# ╔═╡ 31017d7c-d875-442f-b77c-9abc828f42d6
"""
	bidomain(;N=100, dim=1, Δt=1e-4, tₑ=T)

Solves the bidomain problem in `dim` dimensions with `N` grid points in each dimension. Uses Δt as time step until final time tₑ.
"""
function bidomain(;N=100, dim=1, Δt=1e-3, T=30, Plotter=Plots, 
		initial2D=false, Tinit_solve=40, use_csv=false, anisotropic=false)
	xgrid = create_grid(N, dim)

	σᵢ, σₑ = anisotropic ? (σᵢ_anisotropic, σₑ_anisotropic) : (σᵢ_normal, σₑ_normal)
	physics = create_physics(σᵢ,σₑ; anisotropic=anisotropic)
	
	sys = VoronoiFVM.System(
		xgrid,
		physics,
		unknown_storage=:sparse,
	)

	boundaries = (dim == 1 ? 2 : 4)
	enable_species!(sys, species=[1,2,3])
	if !initial2D || true
		boundary_dirichlet!(sys,2,1,0)
	end
	for ispec ∈ [1 3]
		for ibc=1:boundaries
			boundary_neumann!(sys,ispec,ibc,0)
		end
	end

	init = initial_cond(sys, xgrid, Δt, Tinit_solve, 
		dim, initial2D, use_csv, anisotropic)
	U = unknowns(sys)

	SolArray = copy(init)
	tgrid = initial2D ? (Tinit_solve:Δt:T+Tinit_solve) : (0:Δt:T)
	for t ∈ tgrid[2:end]
		solve!(U, init, sys; tstep=Δt)
		init .= U
		SolArray = cat(SolArray, copy(U), dims=3)
	end
	vis = GridVisualizer(resolution=(400,300), dim=dim, Plotter=Plotter)
	return xgrid, tgrid, SolArray, vis, sys
end

# ╔═╡ 8c168042-b17e-434f-b5da-423875d6cc37
species = [L"u", L"u_e", L"v"]

# ╔═╡ 010d60c7-900e-46f4-ae87-2249abd17c16
begin
	dim = 1; N = 1000; Δt = 1e-1;
	xgrid, tgrid, sol, vis = bidomain(dim=dim, N=N, T=T, Δt=Δt);
end;

# ╔═╡ 45c86875-d68a-40d3-b340-d2b4af94e849
function plot_at_t(t,vis,xgrid,sol; title="")
	tₛ = Int16(round(t/Δt))+1
	scalarplot!(vis, xgrid, sol[1,:,tₛ], linestyle=:solid)
	scalarplot!(vis, xgrid, sol[2,:,tₛ], linestyle=:dash)
	plotted_sol = scalarplot!(
		vis, xgrid, sol[3,:,tₛ], linestyle=:dot, legend=:best, show=true, title=title)
	for (i,spec) ∈ zip(2 .* (1:length(species)), species)
		plotted_sol[1][i][:label] = spec
	end
	Plots.plot!(labels=species)
	Plots.plot!(ylims=(-2.5,2.5))
	plotted_sol
end

# ╔═╡ 299c8cf9-4a7e-418c-9bf4-fa8e42da72e6
md"""
#### Solution to 1D problem for t=11
"""

# ╔═╡ d477647c-de1f-4dbf-b982-c1577adcf398
plot_at_t(11,vis,xgrid,sol;title="1D problem at t=11")

# ╔═╡ 6fcf5e8c-c44e-4684-b512-cef6031ad66b
begin
	anim = @animate for t in 0:Δt*5:T
		plot_at_t(t,vis,xgrid,sol)
	end
	gif(anim, "../movies/1D_solution.gif", fps=15)
end;

# ╔═╡ 882e929d-0188-4907-92ef-d9066b92148c
function plot_species_3d(spec)
	Xgrid = xgrid[Coordinates][:]; Tgrid = tgrid; Zgrid = sol[spec,:,:];
	xstep = 1; tstep = 1;
	Tgrid = Tgrid[1:xstep:end]; 
	Xgrid = Xgrid[1:tstep:end]
	Zgrid = Zgrid[1:tstep:end,1:xstep:end]; 
	PyPlot.clf()
	PyPlot.suptitle("Space-Time plot for "*species[spec])
	PyPlot.surf(Xgrid,Tgrid,Zgrid',cmap=:coolwarm) # 3D surface plot
	ax=PyPlot.gca(projection="3d")  # Obtain 3D plot axes
	α = 30; β = 240
	ax.view_init(α,β) # Adjust viewing angles

	PyPlot.xlabel(L"X")
	PyPlot.ylabel(L"T")
	PyPlot.zlabel(L"u")
	figure=PyPlot.gcf()
	figure.set_size_inches(7,7)
	figure
end

# ╔═╡ 61377a01-378d-4ea5-bcfd-b827acfba523
function contour_plot(spec; initial2D=false, anisotropic=false)
	PyPlot.clf()
	Xgrid = xgrid[Coordinates][:]
	Tgrid = tgrid
	Zgrid = sol[spec,:,:]
	PyPlot.suptitle("Space-time contour plot of "*species[spec])
	
	PyPlot.contourf(Xgrid,Tgrid,Zgrid',cmap="hot",levels=100)
	axes=PyPlot.gca()
	axes.set_aspect(2)
	PyPlot.colorbar()

	PyPlot.size(600,600)
	PyPlot.xlabel(L"x")
	PyPlot.ylabel(L"t")
	figure=PyPlot.gcf()
	figure.set_size_inches(7,7)
	extra = initial2D ? ("_inital2D") : (anisotropic ? "_anisotropic" : "") 
	PyPlot.savefig("../img/st_contour_plot_spec_"*string(spec)*extra)
	figure
end

# ╔═╡ 051c235c-ad26-4b29-81c0-6e3e281bb61e
contour_plot(1)

# ╔═╡ 7d2a9f6a-9063-4bc4-bf90-b18dcb30ca9e
contour_plot(2)

# ╔═╡ e7587ca1-7116-4526-8d45-4c76941a98c7
contour_plot(3)

# ╔═╡ 1bda87c6-5e03-4e0a-a029-47d541c3f06a
md"""
### 2D Problem with 1D problem setup
"""

# ╔═╡ a24ebb01-038c-4107-b730-887a86a14fe7
begin
	dim₂=2; N₂ = (100,25); Δt₂ = 1e-1;
	xgrid₂, tgrid₂, sol₂, vis₂ = bidomain(dim=dim₂, N=N₂,T=T, Δt=Δt₂, Plotter=PyPlot);
end;

# ╔═╡ 0edfe024-2786-4570-96f2-f02a083553f6
function contour_2d_at_t(spec, t, Δt, xgrid, sol, title)
	tₛ = Int16(round(t/Δt))+1
	p = scalarplot(
		xgrid,sol[spec,:,tₛ], Plotter=PyPlot, colormap=:hot, 
		title=title, colorlevels=100, isolines=0)
	p.set_size_inches(7,7)
	PyPlot.xlabel(L"x")
	PyPlot.ylabel(L"y")
	p
end

# ╔═╡ 17684f81-f302-4448-b86e-4f0457ad17ac
function contour_subplots(spec, times, xgrid, sol; Δt=1e-1, save=false)
	subplots = [(1,1),(1,2),(2,1),(2,2)]
	p = GridVisualizer(resolution=(700,700),dim=2,Plotter=PyPlot,layout=(2,2),title="Solution with anisotropic conductivity")
	cnt = 0
	for (t,sp) ∈ zip(times,subplots)
		cnt += 1
		tₛ = Int16(round(t/Δt))+1
		scalarplot!(p[sp...],
			xgrid,sol[spec,:,tₛ], Plotter=PyPlot, colormap=:hot, 
			title="t=$t", colorlevels=100, colorbar=false, framepos=cnt)
	end
	p.context[:title] = "Solution with anisotropic conductivity"
	save ? p : reveal(p)
end

# ╔═╡ 4becbe30-d7a1-4949-a494-7f9bba8ed4c9
contour_2d_at_t(2,3,Δt₂,xgrid₂,sol₂,
	"2D problem with 1D problem setup for "*species[2]*" at t="*string(3))

# ╔═╡ b01448fa-c503-4f30-9414-c091bc58b21a
times = [0,10,20,30];

# ╔═╡ 2c2c89f8-5d82-4f43-a679-0b3c5c379f30
md"""
#### Solution to 2D problem with 1D problem setup for $u$:
"""

# ╔═╡ 552af2fb-d17e-4c06-b3f9-856731652f1f
contour_subplots(1,times,xgrid₂,sol₂;Δt=Δt₂)

# ╔═╡ 7dec0320-4406-4360-9e2c-0f0ca2db0afe
md"""
#### Solution to 2D problem with 1D problem setup for $u_e$:
"""

# ╔═╡ 5d46a365-86f1-41a2-b171-600ad4d92549
contour_subplots(2,times,xgrid₂,sol₂;Δt=Δt₂)

# ╔═╡ 9e678ba2-a061-4482-b647-1614494f0b14
md"""
#### Solution to 2D problem with 1D problem setup for $v$:
"""

# ╔═╡ 6ec24d7f-51b4-48f7-8435-5075e4677ee9
contour_subplots(3,times,xgrid₂,sol₂;Δt=Δt₂)

# ╔═╡ 324ecfc7-cdf6-4e2b-9e56-78a9baae0263
md"""
### 2D Problem
"""

# ╔═╡ c1d04237-1696-4c08-9e78-7e7c0d659d0b
begin
	dim₃=2; N₃ = (100,25); Δt₃ = 1e-1;
	xgrid₃, tgrid₃, sol₃, vis₃ = bidomain(
		dim=dim₃, N=N₃, T=T, Δt=Δt₃, Plotter=PyPlot, initial2D=true, Tinit_solve=0)
end;

# ╔═╡ c7106eb9-b54d-4473-8ff4-2b8707260cec
contour_2d_at_t(2,5,Δt₃,xgrid₃,sol₃,
	"2D problem for "*species[2]*" at t="*string(5))

# ╔═╡ c5bbd72e-9e33-496c-8d7a-31b1bf02c00f
md"""
#### Solution to 2D problem for $u$:
"""

# ╔═╡ ac916c5d-0ab5-44d5-8b01-c36c0679ff75
contour_subplots(1,times,xgrid₃,sol₃;Δt=Δt₃)

# ╔═╡ 07a2808d-c680-4118-b9de-23f69369888b
md"""
#### Solution to 2D problem with 1D problem setup for $u_e$:
"""

# ╔═╡ 471473df-ac49-421f-a371-af88e92f14db
contour_subplots(2,times,xgrid₃,sol₃;Δt=Δt₃)

# ╔═╡ a4b916ca-142b-4e26-ad9c-16d5a03acfe9
md"""
#### Solution to 2D problem for $v$:
"""

# ╔═╡ 1dea2df4-2124-4052-bb81-45733cfc5de7
contour_subplots(3,times,xgrid₃,sol₃;Δt=Δt₃)

# ╔═╡ d65e72f0-172d-494f-866f-5b9541035b6f
md"""
### 2D Problem with anisotropic conductivity
"""

# ╔═╡ 36854b9a-39f3-4fca-94da-f9f97cee1663
begin
	dim₄=2; N₄=(100,25); Δt₄=1e-1;
	xgrid₄, tgrid₄, sol₄, vis₄ = bidomain(
		dim=dim₄, N=N₄, T=T, Δt=Δt₄, Plotter=PyPlot, initial2D=true,
		anisotropic=true);
end;

# ╔═╡ 65931b4e-2698-4ddc-990f-3d205bf0c686
contour_2d_at_t(2,21,Δt₄,xgrid₄,sol₄, 
	"2D problem with anisotropic conductivity for "*species[2]*" at t="*string(21))

# ╔═╡ cffaed3f-b4bd-4507-a62f-91831f9729b4
md"""
#### Solution to 2D problem with anisotropic conductivity for $u$:
"""

# ╔═╡ 841f7675-fc4a-43f4-9ab6-1a20caf3315c
contour_subplots(1,times,xgrid₄,sol₄;Δt=Δt₄)

# ╔═╡ 8febdb97-da8f-431b-a0ca-95b857ab138e
md"""
#### Solution to 2D problem with anisotropic conductivity for $u_e$:
"""

# ╔═╡ f41540ac-e2ce-4fbf-9e0c-08a860a9a41d
contour_subplots(2, times, xgrid₄, sol₄; Δt=Δt₄)

# ╔═╡ 22f143da-3dbc-4566-a99d-470fb5e87047
md"""
#### Solution to 2D problem with anisotropic conductivity for $v$:
"""

# ╔═╡ 6c1126d4-33ce-4e10-9a2c-b815b144498b
contour_subplots(3,times,xgrid₄,sol₄;Δt=Δt₄)

# ╔═╡ 81ba9ddf-e84f-463a-9ef7-dc574bf1b940
function save_all_subplots()
	plot_times = [0,10,20,30]
	grids = [xgrid₂,xgrid₃,xgrid₄]
	solutions = [sol₂,sol₃,sol₄]
	Δtₛ = [Δt₂,Δt₃,Δt₄]
	fnames = ["2D_with_1D_setup", "2D", "2D_with_anisotropic"]
	for i ∈ length(grids)
		grid = grids[i]; sol = solutions[i]; Δt = Δtₛ[i]; fname = fnames[i]
		for spec ∈ 1:3
			p = contour_subplots(spec, plot_times, grid, sol; Δt=Δt, save=true)
			GridVisualize.save("../img/"*fname*"_spec"*string(spec)*".png", p)
		end
	end
end

# ╔═╡ 6f9637b0-f2a1-426b-b5cb-a3a7764bcf88
save_all_subplots()

# ╔═╡ 557cc2d0-780a-4f6f-b5c9-6ae9bc31b014
function images_for_gif(spec, folder, xgrid, sol; steps=5)
	folder = "../img/"*folder*"_"*string(spec)*"/"
	mkpath(folder)

	grid_vis = GridVisualizer(resolution=(500,500), dim=2, Plotter=PyPlot)
	for ts=2:steps:size(sol)[3]
		scalarplot!(grid_vis,
			xgrid, sol[spec,:,ts], Plotter=PyPlot, colormap=:hot, colorlevels=100,
			colorbar=false)
		GridVisualize.save(folder*string(ts)*".png", grid_vis)
	end
end

# ╔═╡ 33539ff5-9634-457d-991e-6af44184ce62
for spec=1:3
	images_for_gif(spec, "contour_plot_species", xgrid₂, sol₂; steps=5)
	images_for_gif(spec, "contour_plot_2D_species", xgrid₃, sol₃; steps=5)
	images_for_gif(spec, "contour_plot_2Danisotropic_species", xgrid₄, sol₄; steps=5)
end 	

# ╔═╡ a52838f9-3c48-4653-a7c1-3be6af3109ce
md"""
### Problems
We have not gotten the anisotropic solution or the 2D problem setup to give the same solution as they have in the paper, we never figured out what the problem was. Possible problems includes setting the dirichlet boundary for $u_e$:
```math
\begin{align}
u_e(0,0)=0
\end{align}
```
"""

# ╔═╡ 82e15118-5168-488f-a010-ce1eee2e5a53
md"""
# Performance improvement
To improve performance and stability, as an important first step, we could use more sophisticated time-stepping approaches.

In the Ethier and Borgoult paper, they show a number of schemes and their performance. From this, we can conclude that employing 2nd order methods is most effective. However, we also note that IMEX yields a linear system so each time step, so even if less stable the performance advantage could be very significant when properly implemented.

When it comes to parallelization, the ability to parallelize the solution of the problem at each time step would yield massive performance gains, as it amounts to performing a sequence of sparse matrix-vector multiplications. This is because solving the linear system that results in each iteration of the Newton method can be done by an iterative scheme such as CG, and these schemes also amount to performing a series of sparse matrix-vector multiplications. 
"""

# ╔═╡ a37ba608-ad5d-42ac-a878-bf08aca33c1f
md"""
## Optional topics, Anisotropic conductivity

For the anisotropic conductivity we have different conductivity in each direction. So in the multiplication of $\sigma_i$ and $\sigma_e$ for the fluxes we need to get the contribution to the correct direction. We then need to get the angle between the two edges $x_k$ and $x_l$:
```math
\begin{align}
	\theta = \cos^{-1}\pth{\frac{(x_k-x_l)\cdot e_x}{|x_k-x_l|}}
\end{align}
```
which we then take the corresponding values of the tensor, and we multiply the flux with $\sigma_{i,new}$, $\sigma_{e,new}$ which we get as
```math
\begin{align}
	\sigma_{new} = \sigma_{1,1}\cos(\theta)=\sigma_{2,2}\sin(\theta).
\end{align}
```
"""

# ╔═╡ 6b1e8f04-ae8b-4758-86d9-b947e1620c0d
md"""
# References
"""

# ╔═╡ 510778f6-f368-4bc9-9543-76d0576dbe7a
begin
	references=[
	DOI("10.5281/zenodo.3529808"), # VoronoiFVM
	DOI("10.1016/0025-5564(94)90049-3"), # Bidomain and monodomain introduction
	DOI("10.1137/070680503") #Bidomain problem equations
]
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
ShortCodes = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
CSV = "~0.10.3"
DataFrames = "~1.3.2"
ExtendableGrids = "~0.9.1"
GridVisualize = "~0.5.1"
LaTeXStrings = "~1.3.0"
Plots = "~1.27.0"
PlutoUI = "~0.7.37"
PlutoVista = "~0.8.12"
PyCall = "~1.93.1"
PyPlot = "~2.10.0"
ShortCodes = "~0.3.3"
VoronoiFVM = "~0.16.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "d49f55ff9c7ee06930b0f65b1df2bfa811418475"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.4"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "b15a6bc52594f5e4a3b825858d1089618871bf9d"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.36"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "9310d9495c1eb2e4fa1955dd478660e2ecab1fbb"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.3"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9d3c0c762d4666db9187f363a76b47f7346e673b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.49"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "90b158083179a6ccbce2c7eb1446d5bf9d7ae571"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.7"

[[DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "7eb5d99577e478d23b1ba1faa9f8f6980d34d0a3"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.4"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "a0fcc1bb3c9ceaf07e1d0529c9806ce94be6adf9"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.9"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "Test"]
git-tree-sha1 = "b8abd7c10625b1936f53e8e407c2060a489cc96e"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.9.1"

[[ExtendableSparse]]
deps = ["DocStringExtensions", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "eb3393e4de326349a4b5bccd9b17ed1029a2d0ca"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "0.6.7"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "0dbc5b9683245f905993b51d2814202d75b34f1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.1"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "56956d1e4c1221000b7781104c58c34019792951"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "5d845bccf5d690879f4f5f01c7112e428b1fa543"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "0.5.1"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JSON3]]
deps = ["Dates", "Mmap", "Parsers", "StructTypes", "UUIDs"]
git-tree-sha1 = "8c1f668b24d999fb47baf80436194fdccec65ad2"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.9.4"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "fbd884a02f8bf98fd90c53c1c9d2b21f9f30f42a"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.8.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "4f00cc36fede3c04b8acf9b2e2763decfdcecfa6"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.13"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0886d229caaa09e9f56bcf1991470bd49758a69f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.3"

[[MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "81b44a8cba10ff3cfb564da784bf92e5f834da0e"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.3"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e8185b83b9fc56eb6456200e873ce598ebc7f262"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.7"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "60e9def572717de8345d65a1b913df0fd3903621"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.4"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "9213b4c18b57b7020ee20f33a4ba49eb7bef85e0"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.0"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[PlutoVista]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "GridVisualize", "HypertextLiteral", "UUIDs"]
git-tree-sha1 = "2435d1d3e02db324414f268f30999b5c06a0d10f"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "0.8.12"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "1fc929f47d7c151c839c5fc1375929766fb8edcc"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.93.1"

[[PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "14c1b795b9d764e1784713941e787e1384268103"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.10.0"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "b66df9b4f668b340a6b6b8a7e667a68f586c5561"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.25.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "c086056df381502621dc6b5f1d1a0a1c2d0185e7"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.28.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[ShortCodes]]
deps = ["Base64", "CodecZlib", "HTTP", "JSON3", "Memoize", "UUIDs"]
git-tree-sha1 = "0fcc38215160e0a964e9b0f0c25dcca3b2112ad1"
uuid = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"
version = "0.3.3"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "87efd1676d87706f4079e8e717a7a5f02b6ea1ad"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.2"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "65068e4b4d10f3c31aaae2e6cb92b6c6cedca610"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.6"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "d24a825a95a6d98c385001212dc9020d609f2d4f"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.8.1"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "bfa211c9543f8c062143f2a48e5bcbb226fd790b"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.7"

[[Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "074e08aea1c745664da5c4b266f50b840e528b1c"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.3.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "6dad289fe5fc1d8e907fa855135f85fb03c8fa7a"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.9"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "d60b0c96a16aaa42138d5d38ad386df672cb8bd8"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.16"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c76399a3bbe6f5a88faa33c8f8a65aa631d95013"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.73"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[VoronoiFVM]]
deps = ["DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "IterativeSolvers", "JLD2", "LinearAlgebra", "Parameters", "Printf", "RecursiveArrayTools", "Requires", "SparseArrays", "SparseDiffTools", "StaticArrays", "Statistics", "SuiteSparse", "Symbolics", "Test"]
git-tree-sha1 = "527d3bbded231e626639029d9e6492b2991bd3a0"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.16.2"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═bfef8ade-0b16-4a28-86d7-7e91876519c9
# ╟─c3c70bce-7f23-477b-8f6d-d51b42ed6b2e
# ╟─6c8a63a7-2f5f-4e10-aca5-06d58d0eacd1
# ╟─9653fb91-a348-4c4e-9892-020211969393
# ╟─41ad9bd9-6cea-4b44-a3af-297b21db9f67
# ╟─5c0fa52c-aad8-46c3-a660-4ac85b3904ea
# ╟─b08f1ba2-ff89-4d25-84cf-ef100cde2d91
# ╟─eb03ae91-d8bf-43cd-87b4-9a925b5051f4
# ╟─60a981c1-5e6b-459b-b7b1-dab2b94e1881
# ╟─3ee2e07c-3463-470b-a478-bc7e2464389d
# ╟─2a0bc645-8dec-454e-9eb8-623c7e95a20c
# ╟─2cb6b44b-6236-40e3-a945-4bb8df3cae1e
# ╟─d7aeab82-eab7-44f1-8586-4504287a4dfc
# ╟─042f94d9-04d0-4a03-b3b6-375f8b22953b
# ╟─ac978664-28a5-40bb-8a99-4efbf3c6c6fe
# ╟─63d685c1-dcd1-454f-8e56-f697e7aca8fb
# ╟─b79e3858-2208-4cf9-9834-c9bc1c8e3486
# ╟─3e1c814e-f69d-4839-8ca1-ed03952188a8
# ╠═71943829-1f12-46ce-8347-99a4aebac3d4
# ╠═63372ffa-406b-4eef-a9b2-394326cd35a6
# ╠═d3b4bdc4-f17f-4a83-9ddc-bfdd1ea7f26c
# ╠═d37259c3-f3f7-443c-87fb-497f700e6bac
# ╠═78326e2d-4e56-4b4a-8b4f-6f31662b3f8e
# ╠═b9adb394-034e-4ec5-a539-2a0e850fe1b5
# ╠═7a7229fc-771c-4f87-a705-a29ddc546ac1
# ╟─761e9442-62ab-41da-86e3-cc5ef768dd49
# ╠═47c78a62-d0a1-4120-9303-35a9d7d20ee9
# ╟─dea5cc0a-a112-494a-85ba-0fcf15d6c2e6
# ╟─f4f6357b-3d1b-4730-909c-8a62ae8c284e
# ╠═6cb8de00-62a6-46b0-a367-21a282e4ff89
# ╟─56bec1ef-6b16-4a37-9566-f5fb4e6a2190
# ╠═1fc8a5f0-fa99-44c6-a19a-3334e20eec15
# ╠═31017d7c-d875-442f-b77c-9abc828f42d6
# ╠═6ef697e9-a4b8-446e-bbb8-b577cea0161d
# ╠═d5315aef-978d-447f-85dc-a3becdc33078
# ╠═f19b980c-7d6c-41c6-99af-8475f7aa72db
# ╠═8c168042-b17e-434f-b5da-423875d6cc37
# ╠═010d60c7-900e-46f4-ae87-2249abd17c16
# ╠═45c86875-d68a-40d3-b340-d2b4af94e849
# ╟─299c8cf9-4a7e-418c-9bf4-fa8e42da72e6
# ╠═d477647c-de1f-4dbf-b982-c1577adcf398
# ╠═6fcf5e8c-c44e-4684-b512-cef6031ad66b
# ╠═882e929d-0188-4907-92ef-d9066b92148c
# ╠═61377a01-378d-4ea5-bcfd-b827acfba523
# ╠═051c235c-ad26-4b29-81c0-6e3e281bb61e
# ╠═7d2a9f6a-9063-4bc4-bf90-b18dcb30ca9e
# ╠═e7587ca1-7116-4526-8d45-4c76941a98c7
# ╟─1bda87c6-5e03-4e0a-a029-47d541c3f06a
# ╠═a24ebb01-038c-4107-b730-887a86a14fe7
# ╠═0edfe024-2786-4570-96f2-f02a083553f6
# ╠═17684f81-f302-4448-b86e-4f0457ad17ac
# ╠═4becbe30-d7a1-4949-a494-7f9bba8ed4c9
# ╠═b01448fa-c503-4f30-9414-c091bc58b21a
# ╟─2c2c89f8-5d82-4f43-a679-0b3c5c379f30
# ╠═552af2fb-d17e-4c06-b3f9-856731652f1f
# ╟─7dec0320-4406-4360-9e2c-0f0ca2db0afe
# ╠═5d46a365-86f1-41a2-b171-600ad4d92549
# ╟─9e678ba2-a061-4482-b647-1614494f0b14
# ╠═6ec24d7f-51b4-48f7-8435-5075e4677ee9
# ╟─324ecfc7-cdf6-4e2b-9e56-78a9baae0263
# ╠═c1d04237-1696-4c08-9e78-7e7c0d659d0b
# ╠═c7106eb9-b54d-4473-8ff4-2b8707260cec
# ╟─c5bbd72e-9e33-496c-8d7a-31b1bf02c00f
# ╠═ac916c5d-0ab5-44d5-8b01-c36c0679ff75
# ╟─07a2808d-c680-4118-b9de-23f69369888b
# ╠═471473df-ac49-421f-a371-af88e92f14db
# ╟─a4b916ca-142b-4e26-ad9c-16d5a03acfe9
# ╠═1dea2df4-2124-4052-bb81-45733cfc5de7
# ╟─d65e72f0-172d-494f-866f-5b9541035b6f
# ╠═36854b9a-39f3-4fca-94da-f9f97cee1663
# ╠═65931b4e-2698-4ddc-990f-3d205bf0c686
# ╟─cffaed3f-b4bd-4507-a62f-91831f9729b4
# ╠═841f7675-fc4a-43f4-9ab6-1a20caf3315c
# ╟─8febdb97-da8f-431b-a0ca-95b857ab138e
# ╠═f41540ac-e2ce-4fbf-9e0c-08a860a9a41d
# ╟─22f143da-3dbc-4566-a99d-470fb5e87047
# ╠═6c1126d4-33ce-4e10-9a2c-b815b144498b
# ╠═81ba9ddf-e84f-463a-9ef7-dc574bf1b940
# ╠═6f9637b0-f2a1-426b-b5cb-a3a7764bcf88
# ╠═557cc2d0-780a-4f6f-b5c9-6ae9bc31b014
# ╠═33539ff5-9634-457d-991e-6af44184ce62
# ╟─a52838f9-3c48-4653-a7c1-3be6af3109ce
# ╟─82e15118-5168-488f-a010-ce1eee2e5a53
# ╟─a37ba608-ad5d-42ac-a878-bf08aca33c1f
# ╟─6b1e8f04-ae8b-4758-86d9-b947e1620c0d
# ╟─510778f6-f368-4bc9-9543-76d0576dbe7a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
