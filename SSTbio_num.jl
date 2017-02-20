#hindered settling velocity function
function vhs(X)
  vo/(1.0+(X/Xbar)^q)
end
# Settling batch function
function fb(X)
  X*vhs(X)
end

#Gudonov numerical flux for fb
function Gj(Cl, Cr)
  if Cl < Cr
    min(fb(Cr), fb(Cl))
  else
    max(fb(Cl), fb(Cr))
  end
end

#compression function
function dcomp(X)
  if X < Xc
    0.0
  else
    vhs(X)*ρs*α/(grav*Δρ)
  end
end

#Quadrature integration
using FastGaussQuadrature
const nodes, weights = gausslegendre(10);
function Dcomp(X)
  vals = (X-Xc)/2.0*nodes+(X+Xc)/2.0
  (X-Xc)/2.0*dot(weights, map(dcomp, vals))
end

#Numerical compressive flux
function Jj(Xl, Xr)
  (Dcomp(Xr)-Dcomp(Xl))/Δz
end

#Specific grow rate function
function μj(SNO3, SS)
  μmax*SNO3*SS/((KNO3+SNO3)*(Ks+SS))
end
function PPj(Pl, Pr, Fj)
  if Fj <= 0.0
    Pr
  else
    Pl
  end
end

#Soluble component diffusive fluxes
function Dγj(Sl, Sr)
  ds*(Sr-Sl)/Δz
end

#Update concentrations and P
function update_XP(X, P, SNO3j, Ssj, SNO2j, Xohoj)
  Xn = zeros(X)
  Pn = zeros(P)
  for j = 1:(N+1)
    Xl = j == 1 ? 0.0 : X[j-1]
    Xa = X[j]
    Xr = j==(N+1) ? 0.0 : X[j+1]
    Pa = P[j]
    Pl = j == 1 ? 0.0 : P[j-1]
    Pr = j==(N+1) ? 0.0 : P[j+1]
    Xn[j], Pn[j] = update_conc(Xl, Xa, Xr, Pl, Pa, Pr, SNO3j[j], Ssj[j], Xohoj[j], j)
  end
  return Xn, Pn
end

function update_conc(Xl,Xa,Xr,Pl, Pa, Pr,SNO3,SS,Xoho,j)
  Fl = j == 1 ? 0.0 : Gj(Xl,Xa) - Jj(Xl, Xa)
  Fr = j == (N+1) ? 0.0 : Gj(Xa,Xr) - Jj(Xa, Xr)
  xn = Xa + Δt*(-(Fr-Fl)/Δz+(μj(SNO3,SS)-(1-fp)*b)*Xoho)
  pj = Pa
  if xn > eps()
    PPl = j ==1 ? 0.0 : PPj(Pl,Pa,Fl)
    PPr = j == (N+1) ? 0.0 : PPj(Pa,Pr,Fr)
    pj = 1/xn*(Pa*Xa+Δt*(-(PPr*Fr-PPl*Fl)/Δz+(μj(SNO3,SS)-b)*Xoho))
  end
  return xn, pj
end

# update soluble substrates
function update_solubs(SNO3, Ss, SN2, Xoho)
  SNO3n = zeros(SNO3)
  Ssn = zeros(Ss)
  SN2n = zeros(SN2)
  for j = 1:(N+1)
    SNO3l = j == 1 ? 0.0 : SNO3[j-1]
    SNO3a = SNO3[j]
    SNO3r = j==(N+1) ? 0.0 : SNO3[j+1]
    Ssa = Ss[j]
    Ssl = j == 1 ? 0.0 : Ss[j-1]
    Ssr = j==(N+1) ? 0.0 : Ss[j+1]
    SN2a = SN2j[j]
    SN2l = j == 1 ? 0.0 : SN2[j-1]
    SN2r = j==(N+1) ? 0.0 : SN2[j+1]
    SNO3n[j], Ssn[j], SN2n[j] = update_solubs_j(SNO3l, SNO3a, SNO3r, Ssl, Ssa, Ssr, SN2l, SN2a, SN2r,Xoho[j],j)
  end
  return SNO3n, Ssn, SN2n
end
function update_solubs_j(SNO3l, SNO3a, SNO3r, Ssl, Ssa, Ssr, SN2l, SN2a, SN2r,Xoho,j)
  dγl = j == 1 ? 0.0 : Dγj(SNO3l, SNO3a)
  dγr = j == (N+1) ? 0.0 : Dγj(SNO3a, SNO3r)
  SNO3 = SNO3a + Δt*((dγr-dγl)/Δz-(1.0-Y)/(2.86*Y)*μj(SNO3a,Ssa)*Xoho)
  dγl = j == 1 ? 0.0 : Dγj(Ssl, Ssa)
  dγr = j == (N+1) ? 0.0 : Dγj(Ssa, Ssr)
  Ss = Ssa + Δt*((dγr-dγl)/Δz-(1/Y*μj(SNO3a,Ssa)-(1-fp)*b)*Xoho)
  dγl = j == 1 ? 0.0 : Dγj(SN2l, SN2a)
  dγr = j == (N+1) ? 0.0 : Dγj(SN2a, SN2r)
  SN2 = SN2a + Δt*((dγr-dγl)/Δz+(1.0-Y)/(2.86*Y)*μj(SNO3a,Ssa)*Xoho)
  return SNO3, Ss, SN2
end
