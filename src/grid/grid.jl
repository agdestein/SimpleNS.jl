"""
    Grid{T}()

Nonuniform Cartesian grid with floating point type `T`.
"""
Base.@kwdef struct Grid{T}
    Nx::Int = 10                             # Number of x-volumes
    Ny::Int = 10                             # Number of y-volumes
    xlims::Tuple{T,T} = (0, 1)               # Horizontal limits (left, right)
    ylims::Tuple{T,T} = (0, 1)               # Vertical limits (bottom, top)

    x::Vector{T} = T[]                       # Vector of x-points
    y::Vector{T} = T[]                       # Vector of y-points
    xp::Vector{T} = T[]
    yp::Vector{T} = T[]

    # Number of pressure points in each dimension
    Npx::Int = 0
    Npy::Int = 0

    Nux_in::Int = 0
    Nux_b::Int = 0
    Nux_t::Int = 0
    Nuy_in::Int = 0
    Nuy_b::Int = 0
    Nuy_t::Int = 0

    Nvx_in::Int = 0
    Nvx_b::Int = 0
    Nvx_t::Int = 0
    Nvy_in::Int = 0
    Nvy_b::Int = 0
    Nvy_t::Int = 0

    # Number of points in solution vector
    Nu::Int = 0
    Nv::Int = 0
    NV::Int = 0
    Np::Int = 0

    Ωp::Vector{T} = T[]
    Ω::Vector{T} = T[]

    # For order4
    Ωux::Vector{T} = T[]
    Ωvx::Vector{T} = T[]
    Ωuy::Vector{T} = T[]
    Ωvy::Vector{T} = T[]

    hx::Vector{T} = T[]
    hy::Vector{T} = T[]
    hxi::Vector{T} = T[]
    hyi::Vector{T} = T[]
    hxd::Vector{T} = T[]
    hyd::Vector{T} = T[]

    gxi::Vector{T} = T[]
    gyi::Vector{T} = T[]
    gxd::Vector{T} = T[]
    gyd::Vector{T} = T[]

    Buvy::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)
    Bvux::SparseMatrixCSC{T,Int} = spzeros(T, 0, 0)

    # Separate grids for u, v, and p
    xu::Matrix{T} = zeros(T, 0, 0)
    xv::Matrix{T} = zeros(T, 0, 0)
    yu::Matrix{T} = zeros(T, 0, 0)
    yv::Matrix{T} = zeros(T, 0, 0)
    xpp::Matrix{T} = zeros(T,0, 0)
    ypp::Matrix{T} = zeros(T,0, 0)

    # Ranges
    indu::UnitRange{Int} = 0:0
    indv::UnitRange{Int} = 0:0
    indV::UnitRange{Int} = 0:0
    indp::UnitRange{Int} = 0:0
end

"""
    Grid(Nx, Ny, xlims, ylims; T = eltype(x))

Create nonuniform Cartesian box mesh `x` × `y`.
"""
function Grid(Nx, Ny, xlims, ylims; T = eltype(promote(xlims...)))
    x = LinRange(xlims..., Nx + 1)
    y = LinRange(ylims..., Ny + 1)

    Δx = (xlims[2] - xlims[1]) / Nx
    Δy = (ylims[2] - ylims[1]) / Ny

    # Pressure positions
    xp = (x[1:(end-1)] + x[2:end]) / 2
    yp = (y[1:(end-1)] + y[2:end]) / 2

    # Distance between velocity points
    hx = fill(Δx, Nx)
    hy = fill(Δy, Ny)

    # Distance between pressure points
    gx = zeros(Nx + 1)
    gx[1] = hx[1] / 2
    gx[2:Nx] = (hx[1:(Nx-1)] + hx[2:Nx]) / 2
    gx[Nx+1] = hx[end] / 2

    gy = zeros(Ny + 1)
    gy[1] = hy[1] / 2
    gy[2:Ny] = (hy[1:(Ny-1)] + hy[2:Ny]) / 2
    gy[Ny+1] = hy[end] / 2

    # Number of pressure points
    Npx = Nx
    Npy = Ny
    Np = Npx * Npy

    ## u-volumes
    # x[1]   x[2]   x[3] ....      x[Nx]   x[Nx+1]
    # |      |      |              |       |
    # |      |      |              |       |
    # Dirichlet BC:
    # ULe    u[1]   u[2] ....      u(Nx-1) uRi
    # Periodic BC:
    # u[1]   u[2]   u[3] ....      u[Nx]   u[1]
    # Pressure BC:
    # u[1]   u[2]   u[3] ....      u[Nx]   u[Nx+1]

    # x-dir
    Nux_b = 2               # Boundary points
    Nux_in = Nx             # Inner points
    Nux_t = Nux_in + Nux_b  # Total number

    # Y-dir
    Nuy_b = 2               # Boundary points
    Nuy_in = Ny             # Inner points
    Nuy_t = Nuy_in + Nuy_b  # Total number

    # Total number
    Nu = Nux_in * Nuy_in

    ## v-volumes

    # X-dir
    Nvx_b = 2               # Boundary points
    Nvx_in = Nx             # Inner points
    Nvx_t = Nvx_in + Nvx_b  # Total number

    # Y-dir
    Nvy_b = 2               # Boundary points
    Nvy_in = Ny             # Inner points
    Nvy_t = Nvy_in + Nvy_b  # Total number

    # Total number
    Nv = Nvx_in * Nvy_in

    # Total number of velocity points
    NV = Nu + Nv

    ## Adapt mesh metrics depending on number of volumes

    ## X-direction

    # gxd: differentiation
    gxd = fill(Δx, Nx + 1)

    # hxi: integration and hxd: differentiation
    # Map to find suitable size
    hxi = copy(hx)

    # Restrict Nx+2 to Nux_in+1 points
    xin = x[1:(end-1)]

    hxd = fill(Δx, Nx + 1)

    gxi = copy(hx)

    diagpos = 0

    # Matrix to map from Nvx_t-1 to Nux_in points
    # (used in interpolation, convection_diffusion, viscosity)
    Bvux = spdiagm(Nux_in, Nvx_t - 1, diagpos => ones(Nux_in))

    ## Y-direction

    # gyi: integration and gyd: differentiation
    gyd = fill(Δy, Ny + 1)

    # hyi: integration and hyd: differentiation
    # Map to find suitable size
    hyi = copy(hy)

    # Restrict Ny+2 to Nvy_in+1 points
    yin = y[1:(end-1)]

    hyd = fill(Δy, Ny + 1)

    gyi = copy(hy)

    diagpos = 0

    # Matrix to map from Nuy_t-1 to Nvy_in points
    # (used in interpolation, convection_diffusion)
    Buvy = spdiagm(Nvy_in, Nuy_t - 1, diagpos => ones(Nvy_in))

    ##
    # Volume (area) of pressure control volumes
    Ωp = hyi ⊗ hxi

    # Volume (area) of u control volumes
    Ωu = hyi ⊗ gxi

    # Volume of ux volumes
    Ωux = hyi ⊗ hxd

    # Volume of uy volumes
    Ωuy = gyd ⊗ gxi

    # Volume (area) of v control volumes
    Ωv = gyi ⊗ hxi

    # Volume of vx volumes
    Ωvx = gyi ⊗ gxd

    # Volume of vy volumes
    Ωvy = hyd ⊗ hxi

    Ω = [Ωu; Ωv]

    # Metrics that can be useful for initialization:
    xu = ones(1, Nuy_in) ⊗ xin
    yu = yp ⊗ ones(Nux_in)
    xu = reshape(xu, Nux_in, Nuy_in)
    yu = reshape(yu, Nux_in, Nuy_in)

    xv = ones(1, Nvy_in) ⊗ xp
    yv = yin ⊗ ones(Nvx_in)
    xv = reshape(xv, Nvx_in, Nvy_in)
    yv = reshape(yv, Nvx_in, Nvy_in)

    xpp = ones(Ny) ⊗ xp
    ypp = yp ⊗ ones(Nx)
    xpp = reshape(xpp, Nx, Ny)
    ypp = reshape(ypp, Nx, Ny)

    # Indices of unknowns in velocity vector
    indu = 1:Nu
    indv = Nu .+ (1:Nv)
    indV = 1:NV
    indp = NV .+ (1:Np)

    ## Store quantities in the structure
    params = (;
        Npx,
        Npy,
        Np,
        Nux_in,
        Nux_b,
        Nux_t,
        Nuy_in,
        Nuy_b,
        Nuy_t,
        Nvx_in,
        Nvx_b,
        Nvx_t,
        Nvy_in,
        Nvy_b,
        Nvy_t,
        Nu,
        Nv,
        NV,
        Ωp,
        Ω,
        Ωux,
        Ωvx,
        Ωuy,
        Ωvy,
        hxi,
        hyi,
        hxd,
        hyd,
        gxi,
        gyi,
        gxd,
        gyd,
        Buvy,
        Bvux,
        xu,
        yu,
        xv,
        yv,
        xpp,
        ypp,
        indu,
        indv,
        indV,
        indp,
    )

    Grid{T}(; Nx, Ny, xlims, ylims, x, y, xp, yp, hx, hy, gx, gy, params...)
end
