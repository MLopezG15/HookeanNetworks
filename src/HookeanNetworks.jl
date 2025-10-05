module HookeanNetworks
using LinearAlgebra,Statistics,GLMakie,GeometryBasics

export TriangLattice,ForceCalc,CentroMasa,RecordVideo,ReadState,InnerPolygons,CalcEnergies,RecordVideoWithPolygons

function Distancia(v1,v2) #Distancia entre 2 puntos
    return(norm(v2.-v1)) 
end

function isneigh(R1, R2, v;tol=1e-3) #Checa si hay vecinos 
    for i in v
        if (Distancia(R1+i, R2) < tol)
            return true
        end
    end
    return false
end

function define_network(N1, N2, a1, a2, V, tol=1e-3)

    X = [a1[1]*x+a2[1]*y for x in 0:N1 for y in 0:N2];
    Y = [a1[2]*x+a2[2]*y for x in 0:N1 for y in 0:N2];
    
    edges = [(i, j) 
        for i in eachindex(X)
            for j in 1:length(X) if isneigh([X[i], Y[i]], [X[j], Y[j]], V) ];

    vertices = [X,Y]

    return vertices, edges
end

function acomodar(vertices)
    Ver=zeros(2,length(vertices[1]))
    [Ver[:,i]=[vertices[1][i], vertices[2][i]]   for i in eachindex(vertices[1])]
    return Ver
end

function CentroMasa(Frame)
    Xmasa=mean(Frame[1,:])
    Ymasa=mean(Frame[2,:])
    Masa=[0.0, 0.0]
    Masa[1]=Xmasa;Masa[2]=Ymasa;
    return Masa
end

function CortEnlRand(edges, vertices,k_normal; k_cut=0.0 )
    num_v = size(vertices, 2)
    N = Int(round(sqrt(num_v)))
    idx(u,v) = u + v*N
    targets = Set{Tuple{Int,Int}}()
    interior= Set{Tuple{Int,Int}}()
    for v0 in 0:2:(N-3), u0 in 1:2:(N-2)
        # dos triángulos por super-super celda
        blue = [
            (idx(u0+1,v0),   idx(u0+1,v0+1)),
            (idx(u0+1,v0),   idx(u0,  v0+1)),
            (idx(u0,  v0+1), idx(u0+1,v0+1))
        ]
        red = [
            (idx(u0+1,v0+1), idx(u0+1,v0+2)),
            (idx(u0+2,v0+1), idx(u0+1,v0+2)),
            (idx(u0+1,v0+1), idx(u0+2,v0+1))
        ]
        # elegir 1 enlace al azar de cada triángulo
        for tri in (blue, red)
            for i in tri 
                c,d=i
                push!(interior,(c,d))
            end
            #push!(interior,(c,d))
            a, b = rand(tri)
            a>b && ((a,b)=(b,a))
            push!(targets, (a,b))
        end
    end
    k = fill(k_normal, (size(edges,1),2))
    Kint=[]
    for i in eachindex(k[:,1])
        a, b = Int(edges[i][1]), Int(edges[i][2])
        a>b && ((a,b)=(b,a))
        if (a,b) in interior
            k[i,:]=[k_normal, 25]
            push!(Kint,edges[i])
        end
        if (a,b) in targets
            k[i,:] = [k_cut, 0]
            pop!(Kint)
        end
    end
    return k,Kint
end

function CortEnlOrd(edges, vertices,η,ζ,k_normal; k_cut=0.0, )
    num_v = size(vertices, 2)
    N = Int(round(sqrt(num_v)))
    idx(u,v) = u + v*N
    targets = Set{Tuple{Int,Int}}()
    interior= Set{Tuple{Int,Int}}()
    for v0 in 0:2:(N-3), u0 in 1:2:(N-2)
        # dos triángulos por super-super celda
        blue = [
            (idx(u0+1,v0),   idx(u0+1,v0+1)),
            (idx(u0+1,v0),   idx(u0,  v0+1)),
            (idx(u0,  v0+1), idx(u0+1,v0+1))
        ]
        red = [
            (idx(u0+1,v0+1), idx(u0+1,v0+2)),
            (idx(u0+2,v0+1), idx(u0+1,v0+2)),
            (idx(u0+1,v0+1), idx(u0+2,v0+1))
        ]
        # elegir 1 enlace al azar de cada triángulo
        for tri in (blue, red)
            for i in tri 
                c,d=i
                push!(interior,(c,d))
            end
            a=tri[η][1] # [#Enlace, Triangulo a cortar] 
            b = tri[ζ][2]
            a>b && ((a,b)=(b,a))
            push!(targets, (a,b))
        end
    end
    k = fill(k_normal, (size(edges,1),2))
    Kint=[]
    for i in eachindex(k[:,1])
        a, b = Int(edges[i][1]), Int(edges[i][2])
        a>b && ((a,b)=(b,a))
        if (a,b) in interior
            k[i,:]=[k_normal, (k_normal/2)]
            push!(Kint,edges[i])
        end
        if (a,b) in targets
            k[i,:] = [k_cut, 0]
            pop!(Kint)
        end
    end
    return k,Kint
end

function TriangLattice(N::Int; MethodCut::String="None", η=nothing, ζ=nothing, kval::Float64=1.0)
    a1 = [1, 0]
    a2 = [cos(π/3), sin(π/3)]
    X = [a1[1]*x + a2[1]*y for x in -N÷2:N÷2, y in -N÷2:N÷2]
    Y = [a1[2]*x + a2[2]*y for x in -N÷2:N÷2, y in -N÷2:N÷2]
    V = [a1, a2, [cos(2π/3), sin(2π/3)]]

    B = [(i, j) for i in eachindex(X) for j in 1:length(X) if isneigh([X[i], Y[i]], [X[j], Y[j]], V)]

    vertices, edges = define_network(N, N, a1, a2, V)
    vertices = acomodar(vertices)
    vertices .-= CentroMasa(vertices)

    if MethodCut == "Random"
        Kvec, Kint = CortEnlRand(edges, vertices, kval)
    elseif MethodCut == "Ord"
        if η === nothing || ζ === nothing
            throw(ArgumentError("MethodCut='Ord' requires η and ζ to be provided"))
        end
        Kvec, Kint = CortEnlOrd(edges, vertices, η, ζ, kval)
    elseif MethodCut == "None"
        Kint = []
        Kvec = ones(length(edges), 2) .* kval
    else
        throw(DomainError(MethodCut, "Only 'Random', 'Ord', and 'None' are valid choices"))
    end

    return vertices, edges, Kvec, Kint
end

function wca_force(r,ε, σ)
    r_cut = 2^(1/6) * σ
    r_safe = max(r, 1e-6)
    if r < r_cut
        return 24ε * (2*(σ^12 / (r_safe)^13) - (σ^6 / (r_safe)^7))
    else
        return 0.0
    end
end

function ForceCalc(edges::Vector{Tuple{Int64,Int64}},vertices::Matrix{Float64},Vel::Matrix{Float64},Kvec::Matrix{Float64},t::Float64;m::Float64=1.0,Damp::Bool=true,WCA::Bool=true,GaussPulse::Bool=true,r0::Float64=1.0,γ::Float64=0.2,σF::Float64=0.5,t0::Float64=3.0,A::Float64=1.0,M::Int64=1,ε::Float64=0.1, σ::Float64=0.35,GaussCutOff::Float64=10.0,Thermostat::Bool=true,β=1.0)
    F=zeros(size(vertices))
    r1=[];r2=[];
    for i in eachindex(edges)
        p1=edges[i][1]
        p2=edges[i][2]
        k=Kvec[i,1]
        r1=[vertices[1,p1], vertices[2,p1]]
        r2=[vertices[1,p2], vertices[2,p2]]
        dist=Distancia(r2,r1)
        direccion=(r2.-r1)./dist
        
        mag_elast= -k*(dist-r0)
        F_elas=mag_elast.*direccion 
        fuer_Tot= -1 .*F_elas
        if WCA 
            MagFWCA=wca_force(dist,ε,σ)
            WCAProy= MagFWCA.*(direccion)
            fuer_Tot.-=WCAProy
        end
        F[:,p1].+=fuer_Tot
        F[:,p2].-=fuer_Tot 
    end
    if GaussPulse
        F_img=zeros(size(F))
        DirCm=(vertices[:,M]-CentroMasa(vertices))/norm(vertices[:,M]-CentroMasa(vertices))
        F_img[:,M]=DirCm.*A*exp(-((t-t0)^2)/(2*σF)^2) *(t<GaussCutOff)
        F.+=F_img
    end
    if Damp
        fuer_damp= γ.*Vel
        F.-=fuer_damp 
    end
    if Thermostat
        κ=sqrt(2*m*γ/β)
        F_Lang=κ*randn(size(F))
        F.+=F_Lang
    end
    return F
end

function build_segments(points::Vector{Point2f}, edges::Vector{Tuple{Int, Int}},values)
    seg = [(points[p1], points[p2]) for (p1, p2) in edges]
    colors=values
    return seg,colors
end

function RecordVideo(Sim::Array{Float64,3},Title::String,edges::Vector{Tuple{Int64,Int64}},Kvec::Matrix{Float64};Skips::Int64=10,FR::Int64=50)
    data = [Point2f.(Sim[1, :, t],Sim[2,:,t]) for t in eachindex(Sim[1,1,:])]
    pos = Observable(vec(data[1]));
    disp = Observable(Vector{Vec2f}());  # desplazamientos
    initial = vec(data[1]);  # guarda las posiciones iniciales
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(pos)
    hidedecorations!(ax,ticks=false)
    segments = Observable(Vector{Tuple{Point2f, Point2f}}());
    colors = Observable(Vector{Float32}());
    linesegments!(ax, segments, color=colors,colormap=:vanimo);
    record(fig, "$(Title).gif", 1:Skips:length(Sim[1,1,:])-1; framerate=FR) do t
        current= vec(data[t])
        pos[] = current  # actualiza los puntos
        # Calcula desplazamientos respecto a posición inicial
        next = vec(data[t+1])  # siguiente frame
        disp[] = [40 .*Vec2f(next[j] .- current[j]) for j in eachindex(current)]
        seg,cols=build_segments(current,edges,Kvec[:,2])
        segments[]=seg
        colors[]=cols
    end
end

function angulo(u, v, Frame)
    dx = Frame[1,v] - Frame[1,u]
    dy = Frame[2,v] - Frame[2,u]
    return atan(dy, dx)
end

function build_adj(edges)
    adj = Dict{Int, Vector{Int}}()
    for (u,v) in edges
        push!(get!(adj, u, Int[]), v)
        push!(get!(adj, v, Int[]), u)
    end
    return adj
end

function next_vertex(u, v, adj, Frame; maxdist=1.1)
    ang_uv = angulo(v, u, Frame)
    neigh = adj[v]

    # Filtrar vecinos por distancia ≤ maxdist
    valid_neigh = [w for w in neigh if w != u && Distancia(Frame[:,v,1], Frame[:,w,1]) ≤ maxdist]

    if isempty(valid_neigh)
        return nothing  # no hay siguiente -> camino abierto
    end

    # Calcular ángulos relativos y ordenar
    angs = [(w, mod(angulo(v, w, Frame) - ang_uv, 2π)) for w in valid_neigh]
    sort!(angs, by = x -> x[2])
    return angs[1][1]  # el más cercano en sentido CCW
end

function cleanonecycl(Frame,list,idt)
    picture=Frame[:,:,idt]
    coords=[Point2f.(eachcol(picture[:, poly])) for poly in list]
    return coords
end

function next_vertices(c::Int, u::Int, v::Int, adj, Frame; distCent=sqrt(3)+0.1)
    angcv = mod(angulo(c, v, Frame), 2π)
    neigh = adj[v]
    valid_neigh = [w for w in neigh if w != u &&
        Distancia(Frame[:,c], Frame[:,w]) ≤ distCent]

    if isempty(valid_neigh)
        return Int[]
    end

    angs = [(w, mod(angulo(c, w, Frame) - angcv, 2π)) for w in valid_neigh]
    sort!(angs, by = x -> x[2])
    return [w for (w, _) in angs]
end

function VueltaCentr(c::Int, PPaso, adj, Frame; distCent=sqrt(3)+0.1)
    inicio = PPaso[2]   # primer vértice alrededor del centroide

    function dfs(previo, actual, ciclo)
        candidatos = next_vertices(c, previo, actual, adj, Frame; distCent=distCent)

        for siguiente in candidatos
            nuevo_ciclo = [ciclo; siguiente]

            # si cerramos ciclo válido
            if siguiente == inicio
                return nuevo_ciclo
            end

            # evitar ciclos infinitos
            if length(nuevo_ciclo) > length(keys(adj))
                continue
            end

            res = dfs(actual, siguiente, nuevo_ciclo)
            if !isempty(res)
                return res
            end
        end
        return Int[]  # todos los caminos fallaron → backtrack
    end

    return dfs(inicio, inicio, [inicio])
end

function Hacerpoligonos(Centroides, adj, Frame)
    cycles = Vector{Vector{Int}}()
    for Centr in Centroides
        cycle = VueltaCentr(Centr[1], Centr, adj, Frame)
        if !isempty(cycle) && cycle[1]== cycle[end]
            push!(cycles, cycle)
        end
    end
    return cycles
end

function centroids(Kint,N)
    useful=zeros(Int64,2,length(Kint))
    [useful[:,i]=[Kint[i][1],Kint[i][2]]  for i in eachindex(Kint)]
    return setdiff(Int64(1):Int64(1):Int64((N)^2),useful) ::Vector{Int}
end

function centroidchosen(edges,centroids,N)
    chosen = Tuple{Int, Int}[]
    for c in centroids
        idx = findfirst(t -> c in t && (c+N) in t, edges)
        if !isnothing(idx)
            t = edges[idx]
            reordered = t[1] == c ? t : (t[2], t[1])
            push!(chosen, reordered)
        end
    end
    return chosen
end

function InnerPolygons(Kint::Vector{Any}, edges::Vector{Tuple{Int64,Int64}},Frame::Matrix{Float64})
    N=Int64(sqrt(length(Frame[1,:])))
    adj=build_adj(Kint)
    centroides=centroidchosen(edges,centroids(Kint,N),N)
    cycles = Hacerpoligonos(centroides,adj,Frame)
    return cycles
end

function ReadState(Kint::Vector{Any},Sim::Array{Float64,3},edges::Vector{Tuple{Int64,Int64}})
    cycles=InnerPolygons(Kint,edges,Sim[:,:,1])
    R=zeros(length(cycles),length(Sim[1,1,:]))
    N=Int64(sqrt(length(Sim[1,:,1])))
    for i in eachindex(Sim[1,1,:])
        for (j,c) in enumerate(cycles)
            centr=c[1]-N
            Δϕ=angulo(centr,centr+1-N,Sim[:,:,i])-angulo(centr,centr+1,Sim[:,:,i])
            R[j,i]=Δϕ
        end
    end
    return R
end

function CalcEnergies(edges::Vector{Tuple{Int64,Int64}},Kvec::Matrix{Float64},Frame::Matrix{Float64},VFrame::Matrix{Float64},WCA::Bool=true,r0::Float64=1.0,ε::Float64=0.1, σ::Float64=0.35,m::Float64=1.0)
    K=0.5*m*sum(VFrame[1,:].^2 +VFrame[2,:].^2)
    U=0
    for (i,edge) in enumerate(edges)
        dist=Distancia(Frame[:,edge[1]],Frame[:,edge[2]])
        Δr = dist - r0
        U+= 0.5 * Kvec[i] * Δr^2
        if WCA || dist< 2^(1/6) * σ
            U+= (4ε*((σ/dist)^12 - (σ/dist)^6 )+ε)
        end
    end
    return K,U
end

function RecordVideoWithPolygons(
    Sim::Array{Float64,3},
    Title::String,
    edges::Vector{Tuple{Int,Int}},
    Kvec::Matrix{Float64},
    poly_vertices::Vector{Vector{Int}},
    poly_colors::AbstractMatrix;
    cmap=:vanimo,
    poly_cmap=:plasma,      
    poly_alpha=0.5,         
    Skips::Int=10,
    FR::Int=50
)
    T = size(Sim,3)
    data = [Point2f.(Sim[1,:,t], Sim[2,:,t]) for t in 1:T]

    pos   = Observable(vec(data[1]))
    segs  = Observable(Vector{Tuple{Point2f,Point2f}}())
    cols  = Observable(Vector{Float32}())
    polys = Observable(Vector{GeometryBasics.Polygon}())
    pcols = Observable(Vector{RGBAf}())

    fig = Figure()
    ax = Axis(fig[1,1])
    scatter!(ax, pos; color=:black, markersize=6)
    hidedecorations!(ax, ticks=false)


    polyplot = poly!(ax, polys; color=pcols)
    lineplot = linesegments!(ax, segs; color=cols, colormap=cmap)

    record(fig, "$(Title).gif", 1:Skips:T-1; framerate=FR) do t
        pos[]   = vec(data[t])
        segs[], cols[] = build_segments(pos[], edges, Kvec[:,2])

        # update polygons
        polys[] = [Polygon(Point2f.(Sim[1,p,t], Sim[2,p,t])) for p in poly_vertices]

        Ct = poly_colors[:,t]
        if eltype(Ct) <: Real
            grad = cgrad(poly_cmap)
            nCt = (Ct ) ./ (maximum(Ct)-minimum(Ct))
            pcols[] = [RGBAf(get(grad, x).r, get(grad, x).g, get(grad, x).b, poly_alpha) for x in nCt]
        else
            # assume already colors
            pcols[] = [RGBAf(c.r, c.g, c.b, poly_alpha) for c in Ct]
        end
    end
end

end # module HookeanNetworks
