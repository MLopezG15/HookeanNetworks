# ForceCalc – Cálculo de Fuerzas en Redes de Partículas

`ForceCalc` es una función en Julia para calcular las **fuerzas resultantes** en una red de partículas 2D conectadas por enlaces.  
Incluye interacciones elásticas, amortiguamiento viscoso, repulsión de tipo Weeks–Chandler–Andersen (WCA) y un pulso externo Gaussiano opcional.

---


## Uso

```julia
ForceCalc(
    edges::Vector{Tuple{Int64,Int64}},
    vertices::Matrix{Float64},
    Vel::Matrix{Float64}, #Matriz de Velocidades
    Kvec::Matrix{Float64}, #Lista de vecinos con su valor de constante K, 
    t::Float64;
    Damp::Bool=false,
    WCA::Bool=false,
    GaussPulse::Bool=true,
    r0::Float64=1.0,
    γ::Float64=0.2,
    σF::Float64=0.5,
    t0::Float64=3.0,
    A::Float64=6.0,
    M::Int64=1,
    ε::Float64=0.1,
    σ::Float64=0.35,
    GaussCutOff::Float64=10.0
) -> Matrix{Float64}


