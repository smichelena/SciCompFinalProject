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

"""
    f(u,v)

Returns f(u,v) as in equation (1.4) of Ethier and Bourgalt.
"""
f(u,v) = u - u^3/3 - v

"""
    g(u,v)

Returns g(u,v) as in equation (1.5) of Ethier and Bourgalt.
"""
g(u,v ;γ=0.5, β=1.0) = u + β - γ*v

"""
	fg_zeros()

Returns the approximate solution to f=g=0.
"""
function fg_zeros(;γ=0.5, β=1.0)
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
    u, v = fg_zeros()
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
function u⃗₀(x, y, L)
    return u⃗₀(x, L)
end

"""
    u⃗₀_2D(x,y,L)

Returns vector valued function representing 2D initial conditions for u, uₑ and v 
"""
function u⃗₀_2D(x,y)
    u, v = fg_zeros()
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
	return p
end

"""
	create_physics(σᵢ, σₑ)

Given conductivity tensors σᵢ, σₑ and parameter ε defines physics for problem
"""
function create_physics(σᵢ, σₑ; ε=0.1)
	physics = VoronoiFVM.Physics(
		storage = function(y,u,node)
			y[1] = u[1]
			y[2] = 0
			y[3] = u[3]
		end,
		flux = function(y,u,edge)
			y[1] = -σᵢ*(u[1,2]-u[1,1] + u[2,2]-u[2,1])
			y[2] = σᵢ*(u[1,2]-u[1,1]) + (σᵢ+σₑ)*(u[2,2]-u[2,1])
			y[3] = 0
		end,
		reaction = function(y,u,node)
			y[1] = -f(u[1],u[3])/ε
			y[2] = 0
			y[3] = -ε*g(u[1],u[3])
		end,
	)
	return physics
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
	function bidomain(;N=100, L=70, dim=1, Δt=1e-3, T=T, σᵢ=1.0, σₑ=1.0, Plotter=Plots, dim2_special=false)

Solves bidomain problem.

`N`: number of grid points

`L`: length of domain

`dim`: number of dimensions

`Δt`: time step size for theta scheme

`T`: time to solve for

`σᵢ`: conductivity tensor

`σₑ`: conductivity tensor

`Plotter`: which plotter to use

`dim2_special`: set to true if using 1D data on 2D problem
"""
function bidomain(;N=100, L=70, dim=1, Δt=1e-3, T=T, σᵢ=1.0, σₑ=1.0, Plotter=Plots, dim2_special=false)

	xgrid = create_grid(N, L, dim)

	physics = create_physics(σᵢ, σₑ)

	sys = VoronoiFVM.System(
		xgrid,
		physics,
		unknown_storage=:sparse,
	)

    impose_neumann_boundary(sys, dim)

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
	return plotted_sol
end

"""
	plot_species_3d(spec, sol, label)

Plots solucion to 1D bidomain problem as a 3D plot
"""
function plot_species_3d(spec, sol, label, xgrid, tgrid)
	Xgrid = xgrid[Coordinates][:]; Tgrid = tgrid; Zgrid = sol[spec,:,:];
	xstep = 1; tstep = 1;
	Tgrid = Tgrid[1:xstep:end]; 
	Xgrid = Xgrid[1:tstep:end]
	Zgrid = Zgrid[1:tstep:end,1:xstep:end]; 
	PyPlot.clf()
	PyPlot.suptitle("Space-Time plot for "*label[spec])
	PyPlot.surf(Xgrid,Tgrid,Zgrid',cmap=:coolwarm) # 3D surface plot
	ax=PyPlot.gca(projection="3d")  # Obtain 3D plot axes
	α = 30; β = 240
	ax.view_init(α,β) # Adjust viewing angles
	
	PyPlot.xlabel(L"X")
	PyPlot.ylabel(L"T")
	PyPlot.zlabel(L"u")
	figure=PyPlot.gcf()
	figure.set_size_inches(7,7)
	return figure
end

"""
	contour_plot(spec, sol, label)

Makes a contour plot of species spec
"""
function contour_plot(spec, sol, label, xgrid, tgrid)
	PyPlot.clf()
	Xgrid = xgrid[Coordinates][:]
	Tgrid = tgrid
	Zgrid = sol[spec,:,:]
	PyPlot.suptitle("Space-time contour plot of "*label[spec])
	
	PyPlot.contourf(Xgrid,Tgrid,Zgrid',cmap="hot",levels=100)
	axes=PyPlot.gca()
	axes.set_aspect(2)
	PyPlot.colorbar()

	PyPlot.size(600,600)
	PyPlot.xlabel(L"x")
	PyPlot.ylabel(L"t")
	figure=PyPlot.gcf()
	figure.set_size_inches(7,7)
	figure
end

"""
	contour_2d_at_t(spec, t, xgrid, sol, label)

plots contour of 2D priblem at time t
"""
function contour_2d_at_t(spec, t, xgrid, sol, label)
	tₛ = Int16(round(t/Δt₂))+1
	p = scalarplot(
		xgrid,sol[spec,:,tₛ], Plotter=PyPlot, colormap=:hot, 
		title="2D problem with 1D problem setup for "*label[spec]*
		" at t="*string(t),
		colorlevels=100, isolines=0)
	p.set_size_inches(7,7)
	PyPlot.xlabel(L"x")
	PyPlot.ylabel(L"y")
	p
end

function main()
	# β = 1.0
	# γ = 0.5
	# ε = 0.1
	# σᵢ = 1.0
	# σₑ = 1.0
    # L = 70
	# T = 30

	#this is universal so I am leaving it here
	species = [L"u", L"u_e", L"v"]

	#plot intial conditions
    #plot_initial_conditions_1D(70)

	#solve problem in 1D
    xgrid, tgrid, sol, vis = bidomain(dim=1, N=1000, Δt=1e-1);
    #plot_at_t(11, 1e-1, vis, xgrid, sol, species)

	#plot_species_3d(1, sol, species, xgrid, tgrid)

	dim₂=2; N₂ = (100,25); Δt₂ = 1e-1;
	xgrid₂, tgrid₂, sol₂, vis₂ = bidomain(dim=dim₂, N=N₂, Δt=Δt₂, Plotter=PyPlot);
	contour_2d_at_t(2,3,xgrid₂,sol₂, species)
    
end