# HookeanNetworks -Paquete para construcción y simulación de redes elasticas
Actualmente consta de las principales funciones: [TriangLattice](#TriangLattice), [ForceCalc](#ForceCalc), [CentroMasa](#CentroMasa), [RecordVideo](#RecordVideo), [InnerPolygons](#InnerPolygons), [ReadState](#ReadState), [CalcEnergies](#CalcEnergies), [RecordVideoWithPolygons](#RecordVideoWithPolygons)

## TriangLattice
`TriangLattice` es la función que permite obtener un arreglo de enlaces y vertices en un lattice triangular

### Argumentos
- `N::Int64` Número de veces a aplicar el vector de la base, generando un arreglo (N+1)x(N+1)

### Keywords Opcionales
- `MethodCut::String="None"` Método para cortar enlaces internos de la supercelda `MethodCut="Random"` elige cortar enlaces de manera aleatoria. `MethodCut="Ord"` genera cortes periodicos basado en la elección de enlace a cortar.
- `η=nothing` Indice del enlace a cortar para el primer triangulo.
- `ζ=nothing` Indice del enlace a cortar para el segundo triangulo.
- `kval::Float64=1.0` Valor de la constante de resorte, default es 1.

### Uso 

```julia
TriangLattice(N::Int64; MethodCut::String="None", η=nothing, ζ=nothing, kval::Float64=1.0)
```

## ForceCalc

`ForceCalc` es una función en Julia para calcular las **fuerzas resultantes** en una red de partículas 2D conectadas por enlaces.  
Incluye interacciones elásticas, amortiguamiento viscoso, repulsión de tipo Weeks–Chandler–Andersen (WCA) y un pulso externo Gaussiano opcional.

---
### Argumentos
- `edges::Vector{Tuple{Int64,Int64}}` Aristas del arreglo.
- `vertices::Matrix{Float64}` Posiciones de los vertices, de tamaño DxN con D= Número de dimensiones, N= Número de vertices.
- `Vel::Matrix{Float64}` Velocidades en el instante de tiempo que se calcula la fuerza.
- `Kvec::Matrix{Float64}` Matriz Mx2 M:Número de aristas del arreglo. La primera columna da la contante k del resorte, la segunda hace distinción entre los enlaces internos formados por superceldas de 4 triangulos.
- `t::Float64` Valor del tiempo.

### Keywords Opcionales
- `r0::Float64=1.0` Distancia de reposo para los enlaces de la red.
- `Damp::Bool=true` Activa el amortiguamientos en los resortes.
- `γ::Float64=0.2` Viscosidad del amortiguamiento.
- `WCA::Bool=true` Activa el potencial Weeks-Chandler-Andersen para repulsión.
- `ε::Float64=0.1` Parámetro de interacción WCA.
- `σ::Float64=0.35` Parámetro de distancia WCA.
- `GaussPulse::Bool=true` Activa una fuerza externa que actúa como una gaussiana sobre 1 particula.
- `M::Int64=1` Indice de la particula a aplicar la fuerza externa.
- `σF::Float64=0.5` Ancho del pulso gaussiano.
- `t0::Float64=3.0` Desplazamiento del pulso gaussiano.
- `A::Float64=1.0` Amplitud del pulso gaussiano.
- `GaussCutOff::Float64=10.0` Tiempo máximo en el que se aplica el pulso gaussiano.
- `Thermostat::Bool=true` Activa el termostáto de Langevin.
- `β` Temperatura reducida.
- `m` Masa de los vertices del sistema.

### Uso

```julia
ForceCalc(
    edges::Vector{Tuple{Int64,Int64}},vertices::Matrix{Float64},Vel::Matrix{Float64},Kvec::Matrix{Float64},t::Float64;
    Damp::Bool=false,WCA::Bool=false,GaussPulse::Bool=true,r0::Float64=1.0,γ::Float64=0.2,σF::Float64=0.5,t0::Float64=3.0,
    A::Float64=6.0,M::Int64=1,ε::Float64=0.1,σ::Float64=0.35, GaussCutOff::Float64=10.0,Thermostat::Bool=true,β=1.0) -> Matrix{Float64}
```

## CentroMasa

`CentroMasa`es una función que calcula y devuelve la posición del centro de masa del sistema para un instante de tiempo.

### Argumentos
`Frame::Matrix{Float64}` Matriz de posiciones de los vertices del sistema.

### Uso
```julia
CentroMasa(Frame::Matrix{Float64})-> Vector{Float64}
```

## RecordVideo
`RecordVideo` es una función que automáticamente graba en un archivo `.gif` la dinámica de la red. 

### Argumentos
`Sim::Array{Float64,3}` Arreglo de 3 dimensiones (N,D,T) con la posición de las N particulas en D dimensiones a lo largo de un tiempo T.
`Title::String` Titulo del archivo `.gif`.

### Keywords Opcionales
`Skips::Int64` Especialmente útil para simulaciones largas, intervalos de tiempo sobre los que se grafican.
`FR::Int64` FrameRate de la simulación. 

## InnerPolygons
`InnerPolygons` es una función que obtiene las trayectorias cerradas formadas por los enlaces internos de la red.
### Argumentos
- `Kint::Vector{Any}` Vector de Tuple con los indices de los vertices que forman enlaces internos del material.
-  `edges::Vector{Tuple{Int64,Int64}}` Vector con los enlaces de toda la red.
-  `Frame::Matrix{Float64}` Matriz de posición de todos vertices en un tiempo específico.

### Uso
```julia
InnerPolygons(Kint::Vector{Any}, edges::Vector{Tuple{Int64,Int64}},Frame::Matrix{Float64})-> Vector{Vector{Float64}}
```

## ReadState
`ReadState` obtiene a lo largo del tiempo el ángulo de apertura de los modos suaves.
### Argumentos
- `Kint::Vector{Any}` Vector de Tuple con los indices de los vertices que forman enlaces internos del material.
- `Sim::Array{Float64,3}` Arreglo de 3 dimensiones (N,D,T) con la posición de las N particulas en D dimensiones a lo largo de un tiempo T.
- `edges::Vector{Tuple{Int64,Int64}}` Vector con los enlaces de toda la red.
### Uso
```julia
ReadState(Kint::Vector{Any},Sim::Array{Float64,3},edges::Vector{Tuple{Int64,Int64}})-> Matrix{Float64}
```

## CalcEnergies
`CalcEnergies` calcula por separado la energía cinética y energía potencial para un instante de tiempo.
### Argumentos
- `edges::Vector{Tuple{Int64,Int64}}` Vector con los enlaces de toda la red.
- `Kvec::Matrix{Float64}` Matriz Mx2 M:Número de aristas del arreglo. La primera columna da la contante k del resorte, la segunda hace distinción entre los enlaces internos formados por superceldas de 4 triangulos.
- `Frame::Matrix{Float64}`  Matriz de posiciones de todos vertices en un tiempo específico.
- `VFrame::Matrix{Float64}`Matriz de velocidades de todos vertices en un tiempo específico.
### Keywords Opcionales
- `WCA::Bool=true` Considera el potencial WCA para el cálculo de energía total.
- `r0::Float64=1.0` Distancia de equilibrio de los enlaces.
- `ε::Float64=0.1` Parámetro de interacción WCA.
- `σ::Float64=0.35` Parámetro de distancia WCA.
- `m::Float64=1.0` Masa de los vertices del sistema.
### Uso
```julia
CalcEnergies(edges::Vector{Tuple{Int64,Int64}},Kvec::Matrix{Float64},Frame::Matrix{Float64},VFrame::Matrix{Float64},WCA::Bool=true,r0::Float64=1.0,ε::Float64=0.1, σ::Float64=0.35,m::Float64=1.0) -> Float64, Float64
```

## RecordVideoWithPolygons

### Argumentos 
- `Sim::Array{Float64,3}` Arreglo de 3 dimensiones (N,D,T) con la posición de las N particulas en D dimensiones a lo largo de un tiempo T.
- `Title::String`
- `edges::Vector{Tuple{Int,Int}}` Vector con los enlaces de toda la red.
- `Kvec::Matrix{Float64}` Matriz Mx2 M:Número de aristas del arreglo. La primera columna da la contante k del resorte, la segunda hace distinción entre los enlaces internos formados por superceldas de 4 triangulos.
- `poly_vertices::Vector{Vector{Int}}` Vector que contiene la lista de vertices que forman poligonos cerrados.
- `poly_colors::AbstractMatrix` Matriz MxT de valores RGBA para M poligonos a lo largo de un tiempo T.


### Keywords Opcionales
- `cmap=:vanimo` Colormap de los enlaces.
- `poly_cmap=:plasma` Colormap de los polígonos.
- `poly_alpha=0.5` Transparencia de los polígonos.
- `Skips::Int64=10` Especialmente útil para simulaciones largas, intervalos de tiempo sobre los que se grafican.
- `FR::Int64=50` FrameRate de la simulación. 

### Uso
```julia
RecordVideoWithPolygons(Sim::Array{Float64,3},Title::String,edges::Vector{Tuple{Int,Int}},Kvec::Matrix{Float64},poly_vertices::Vector{Vector{Int}},    poly_colors::AbstractMatrix;cmap=:vanimo,poly_cmap=:plasma,poly_alpha=0.5,Skips::Int=10,FR::Int=50)
```


