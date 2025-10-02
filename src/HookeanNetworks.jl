module HookeanNetworks
using LinearAlgebra,Statistics,GLMakie

export TriangLattice,ForceCalc,CentroMasa,RecordVideo

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
    #CCM=UnCentroMasa(vertices)
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

function RecordVideo(Sim::Array{Float64,3},Title::String,Skips::Int64=10,FR::Int64=50)
    data = [Point2f.(Sim[1, :, t],Sim[2,:,t]) for t in eachindex(T)]
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
    record(fig, "$(Title).gif", 1:Skips:length(T)-1; framerate=FR) do t
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

end # module HookeanNetworks
