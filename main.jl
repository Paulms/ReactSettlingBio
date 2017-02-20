# Parameter values employed for the simulation
const vo = 0.00176 #m/s
const Xbar = 3.87 #Kg/m³
const q = 3.58 #-
const fp = 0.2 #-
const Xc = 5 #Kg/m³
const α = 0.2 #m²/2²
const ρs = 1050 #Kg/m³
const Δρ = 52 #Kg/m³
const grav = 9.81 #m/s²
const ds = 0.000001 #m²/s
const Y = 0.67 #-
const μmax = 5.56e-9 #1/s
const b = 6.94e-6 #1/s
const Ks = 0.02 #Kg/m³
const KNO3 = 5e-4 #Kg/m³

#Dimensiones del tanque
const B = 1

#Discretización
const N = 100
const Δz = B/N
const nmax= 100000    # TODO

#Initial conditions
const Sso = 9e-4 #Kg/m³
const SNO3o = 6e-3 #Kg/m³
const Xo = 3.5 #Kg/m³
const po = 5/7

#Load numerical method functions
include("SSTbio_num.jl")

function initialize(M)
  #Variables
  Xj = zeros(N+1, M)

  #Step 2: Initialize by Averaging
  for j = 1:(N+1)
    Xj[j, 1] = Xo#*((j-0.5)*Δz)
  end
  Pj = ones(N+1,M)*po
  Xohoj = po*Xj
  SNO3j = SNO3o*ones(N+1,M)
  Ssj = Sso*ones(N+1,M)
  SN2j = zeros(N+1,M)
  return Xj, Pj, Xohoj, SNO3j, Ssj, SN2j
end

#Step 3: Time iteration
function simular(Tmax, Xjm, Pjm, Xohojm, SNO3jm, Ssjm, SN2jm)
  #Iteramos sobre el tiempo
  Xj = Xjm[:,1]
  Pj = Pjm[:,1]
  Xohoj = Xohojm[:,1]
  SNO3j = SNO3jm[:,1]
  Ssj = Ssjm[:,1]
  SN2j = SN2jm[:,1]
  tiempo = 0
  M = 1
  while tiempo <= Tmax
    tiempo = tiempo + Δt
    if tiempo%60 < Δt
      println("Tiempo: $tiempo seg.")
      M = M + 1
      Xjm[:,M] = Xj
      Pjm[:,M] = Pj
      Xohojm[:,M] = Xohoj
      SNO3jm[:,M] = SNO3j
      Ssjm[:,M] = Ssj
      SN2jm[:,M] = SN2j
    end
    #Update Layers
    Xj, Pj = update_XP(Xj, Pj, SNO3j, Ssj, SN2j,Xohoj)
    SNO3j, Ssj, SN2j = update_solubs(SNO3j, Ssj, SN2j, Xohoj)
    Xohoj = Pj.*Xj
  end
end

#simulate
const Tmax = 3600*2#3600 #h -> s
const Δt = Tmax/nmax
const MM = Int(ceil(Tmax/60)+1)  #Save each minute
println("Starting...")
Xj, Pj, Xohoj, SNO3j, Ssj, SN2j = initialize(MM)
println("Processing...")
@time simular(Tmax, Xj, Pj, Xohoj, SNO3j, Ssj, SN2j)
println("Finished...")

#Plots
using PyPlot
fig = figure(figsize=(8,6))
ax = fig[:gca](projection="3d")
Z = (0:N)*B/N
grid_a = [i for i in Z, j in 1:MM]
grid_b = [j for i in Z, j in 1:MM]
ax[:plot_surface](grid_a, grid_b,Xj, rstride=2, cstride=2, cmap=ColorMap("jet"), alpha=0.7, linewidth=0.25)
xlabel("z [m]")
ylabel("t [min]")
zlabel("X(z,t) [Kg/m³]")
fig[:show]()

#Profile
#Profile.init(delay=0.01)
#Profile.clear()
#@profile simular(Tmax, Xj, Pj, Xohoj, SNO3j, Ssj, SN2j)
#using ProfileView
#ProfileView.view()
