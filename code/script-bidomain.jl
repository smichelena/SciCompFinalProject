#includes
using ShortCodes
using ExtendableGrids
using VoronoiFVM
using PlutoVista
using Plots
using GridVisualize
using PlutoUI
using LaTeXStrings
using PyPlot
using PyCall

#package setup
GridVisualize.default_plotter!(Plots);
pyplot();

begin
	β = 1.0
	γ = 0.5
	ε = 0.1
	σᵢ = 1.0
	σₑ = 1.0
    L = 70
	T = 30
end;

"""
    f(u,v)

Returns f(u,v) as in equation (1.4) of Ethier and Bourgalt.
"""
f(u,v) = u - u^3/3 - v

"""
    g(u,v)

Returns g(u,v) as in equation (1.5) of Ethier and Bourgalt.
"""
g(u,v) = u + β - γ*v

"""
	fg_zeros()

Returns the approximate solution to f=g=0.
"""
function fg_zeros(β,γ)
	p = 3(1-γ)/γ; q = 3β/γ;
	h(u) = u^3 + p*u + q; ḣ(u) = 3u^2 + p
	uₙ = 1; hₙ = 1
	cnt = 0
	while abs(hₙ) > 1e-15
		uₙ = uₙ - h(uₙ)/ḣ(uₙ)
		hₙ = h(uₙ)
	end
	vₙ = uₙ - uₙ^3/3
	return uₙ, vₙ
end

"""
	create_grid(N, dim=1)

Creates a simplex grid in dimension `dim` with `N` grid points in one
dimension with length `L`.
"""
function create_grid(N, L, dim=1)
	if dim == 2
		X = LinRange(0,L,N[1])
		Y = LinRange(0,L,N[2])
		return simplexgrid(X,Y)
	else
		X = LinRange(0,L,N)
		return simplexgrid(X)
	end
end

"""
    u⃗₀(x,L)

Returns vector valued function representing 1D initial conditions for u, uₑ an v
"""
function u⃗₀(x, L)
    u, v = fg_zeros(β, γ)
    if 0 ≤ x ≤ L/20
        u = 2
    end
    uₑ = 0
    return [u, uₑ, v]
end

"""
    u⃗₀(x,y,L)

Returns vector valued function representing 2D initial conditions for u, uₑ and v using 1D data
"""
function u⃗₀(x,y, L)
    return u⃗₀(x, L)
end

"""
    u⃗₀_2D(x,y,L)

Returns vector valued function representing 2D initial conditions for u, uₑ and v 
"""
function u⃗₀_2D(x,y)
    u, v = fg_zeros(β, γ)
    if 0 ≤ x ≤ 3.5 && 0 ≤ y ≤ 70
        u = 2
    end
    if 31 ≤ x ≤ 39 && 0 ≤ y ≤ 35
        v = 2
    end
    uₑ = 0
    return [u, uₑ, v]
end

"""
    plot_initial_conditions_1D(L)

Plots initial conditions for 1D bidomain problem
"""
function plot_initial_conditions_1D(L)
    X=0:0.1:L
	initial = u⃗₀.(X, L)
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

"""
    Define physics for 1D problem
"""
physics = VoronoiFVM.Physics(
	storage = function(y,u,node)
		y[1] = u[1]
		y[2] = 0
		y[3] = u[3]
	end,
	flux = function(y,u,edge)
		nspecies=3
		y[1] = -σᵢ*(u[1,2]-u[1,1] + u[2,2] - u[2,1])
		y[2] = σᵢ*(u[1,2]-u[1,1]) + (σᵢ+σₑ)*(u[2,2]-u[2,1])
		y[3] = 0
	end,
	reaction = function(y,u,node)
		node.index
		y[1] = -f(u[1],u[3])/ε
		y[2] = 0
		y[3] = -ε*g(u[1],u[3])
	end,
)

"""
    bc(y, u, node)

imposes boundary conditions (?)
"""
function bc(y,u,node)
    boundary_dirichlet!(y,u,node; region=1, value=0, species=1)
    for i=1:2
        for j=1:2
            boundary_neumann!(y,u,node; region=i, value=0, species=j)
        end
    end
end

"""
    impose_neumann_boundary(sys, dim)

imposes neumann boundary conditions
"""
function impose_neumann_boundary(sys, dim)
    boundaries = (dim == 1 ? 2 : 4)
    enable_species!(sys, species=[1,2,3])
    boundary_dirichlet!(sys,2,1,0)
    for ispec ∈ [1 3]
        for ibc=1:boundaries
            boundary_neumann!(sys,ispec,ibc,0)
        end
    end
end

"""
	solve_in_time(U, sys, init, inival, t₀, Δt, T; storing=storing)

Solves parabolic PDE until time T
"""
function solve_in_time(U, sys, init, inival, t₀, Δt, T; storing=storing)
	init .= [tuple[k] for tuple in inival, k in 1:3]'
    if !storing
		for t ∈ t₀:Δt:T
			solve!(U, init, sys; tstep=Δt)
			init .= U
		end
        return init
    end
    SolArray = copy(init)
    for t ∈ (t₀ + Δt):Δt:T
		solve!(U, init, sys; tstep=Δt)
		init .= U
		SolArray = cat(SolArray, copy(U), dims=3)
	end
    return SolArray
end

"""
	bidomain(;N=100, dim=1, Δt=1e-4, tₑ=T)

Solves the bidomain problem in `dim` dimensions with `N` grid points in each dimension. 
Uses Δt as time step until final time tₑ.
"""
function bidomain(;N=100, L=70, dim=1, Δt=1e-3, T=T, Plotter=Plots, dim2_special=false)
	xgrid = create_grid(N, L, dim)


	sys = VoronoiFVM.System(
		xgrid,
		physics,
		unknown_storage=:sparse,
	)

    impose_neumann_boundary(sys, dim)

    #this can be factored out
	init = unknowns(sys)
	U = unknowns(sys)
	if dim==2 && dim2_special
		inival = map(u⃗₀_2D, xgrid)
		init .= solve_in_time(U,sys, init, inival, 0, Δt, T; storing=false)
	else
		inival = map(x -> u⃗₀(x,L), xgrid)
		init .= [tuple[k] for tuple in inival, k in 1:3]'
	end

	U = unknowns(sys)
	t₀ = dim2_special ? T : 0
    T = dim2_special ? 2*T : T
	SolArray = solve_in_time(U, sys, init, inival, t₀, Δt, T; storing=true)

	vis = GridVisualizer(resolution=(400,300), dim=dim, Plotter=Plotter)

	return xgrid, t₀:Δt:T, SolArray, vis, sys
end

"""
    plot_at_t(t, vis, xgrid, sol)

Plots solution at time t
"""
function plot_at_t(t, Δt, vis, xgrid, sol, species)
	tₛ = Int16(round(t/Δt))+1
	scalarplot!(vis, xgrid, sol[1,:,tₛ], linestyle=:solid)
	scalarplot!(vis, xgrid, sol[2,:,tₛ], linestyle=:dash)
	plotted_sol = scalarplot!(
		vis, xgrid, sol[3,:,tₛ], linestyle=:dot, legend=:best, show=true)
	for (i,spec) ∈ zip(2 .* (1:length(species)), species)
		plotted_sol[1][i][:label] = spec
	end
	Plots.plot!(labels=species)
	Plots.plot!(ylims=(-2.5,2.5))
	plotted_sol
end

function main()
    plot_initial_conditions_1D(70)

	species = [L"u", L"u_e", L"v"]
    dim = 1; N = 1000; Δt = 1e-1;
    xgrid, tgrid, sol, vis = bidomain(dim=dim, N=N, Δt=Δt);
    plot_at_t(11, Δt, vis, xgrid, sol, species)
    
end