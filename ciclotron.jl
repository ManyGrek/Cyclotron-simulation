using Plots: print

## Vacas Esféricas©

using Plots
using LinearAlgebra

## Parámetros iniciales

const q = 1.602176487e-19   # Carga del protón [C] 
m = 1.6726e-27              # Masa del protón en [kg]

radio = 5e-2;           # Radio del ciclotrón [m]
const d = 0.5e-2;       # Separación entre las placas [m]
const Bz = 1;           # Campo magnético [T]
V = 5000;               # Diferencia de potencial máxima [Volts]
const E = V/d; # Magnitud del campo eléctrico

ω = q*[0,0,Bz]/m; #Frecuencia del ciclotrón requerida
f = 
print("Frecuencia angular ω del ciclotrón: ",norm(ω))


dt = 1e-10            # Paso de la Simulación [s]
tf = (2π/norm(ω))*100 # Tiempo final de la Simulación [s]
t = 0:dt:tf           # Vector de tiempo

# Condiciones iniciales del protón
r0 = [0 0 0]; 
v0 = [0 0 0]; 

## Vectorización de las variables

r = zeros(length(t),3); # Vector de posición
v = zeros(length(t),3); # Vector de velocidad
a = zeros(length(t),3); # Vector de aceleración

r[1,:] = r0; # Aplicando condiciones iniciales a la trayectoria
v[1,:] = v0; # Aplicando condiciones iniciales a la velocidades

## Función de la Fuerza en el tiempo
function F(r,v)
    if abs(r[1]) ≤ d/2
        return q*([sign(v[1])*E,0,0] + cross(v,[0,0,Bz]));
    else
        return q*(cross(v,[0,0,Bz]));
    end
end

## Algoritmo Euler-Cromer
#Primer paso
Fuerza = q*([E,0,0]);

r[2,:] = r[1,:] + dt*v[1,:];
v[2,:] = v[1,:] + dt*Fuerza/m;

for i = 3:length(t)
    r[i,:] = r[i-1,:] + dt*v[i-1,:];
    v[i,:] = v[i-1,:] + dt*F(r[i-1,:],v[i-1,:])/m;
    if (r[i,1] - sign(r[i,1])*d/2)^2+r[i,2]^2 > radio^2
        global paso_final = i;
        break
    end
end

# Eliminar datos inutiles

t=t[1:paso_final,:];
r=r[1:paso_final,:];
v=v[1:paso_final,:];

## Gráfica de la Simulación

plot(r[:,1],r[:,2];
xlabel="x[m]",
ylabel="y[m]",
aspect_ratio=:equal,
legend=false,
xlims=(-0.06,0.06),
ylims=(-0.06,0.06),
title="Trayectoria de un protón en un ciclotrón")

## Gráficas de las Ds #

TH=collect(-π/2:π/128:π/2)
XD=radio*[cos(θ) for θ in TH];
YD=radio*[sin(θ) for θ in TH];

p1=[d/2,d/2];
p2=[-radio,radio];

plot!(XD+d/2 * ones(length(TH),1),YD)
plot!(-XD-d/2 * ones(length(TH),1),-YD)
plot!(p1,p2)
plot!(-p1,-p2)


## Gráfica de le energía cinética K

# Función que arroja la magnitud de cada vector en los
# renglones de una matriz
magnitud(A) = [norm(A[i,:]) for i=1:size(A,1)] 
v_magnitud = magnitud(v)

K = (0.5*m.*v_magnitud.^2)*6.242e18; #Energía cinética en eV

plot(t,K;
xlabel="t[s]",
ylabel="E[eV]",
legend=false,
title="Energía cinética a lo largo del tiempo")


## Tiempo final
tf = dt*paso_final;
print("\nTiempo del protón en el borde exterior del ciclotrón: ",tf," s")

## Energía cinética en el borde exterior
print("\nEnergía cinética en el borde exterior: ",K[paso_final]," eV")

## Velocidad final del protón
print("\nVelocidad final del protón: ",v_magnitud[paso_final]," m/s")