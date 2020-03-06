using DifferentialEquations
using Random
using PyPlot
using Printf
using Statistics
using DataFrames
using IterableTables
using LinearAlgebra
using DelimitedFiles
using CSV
using LinearAlgebra


function regular(u, n, A)
    G=zeros((n,n))
    for i=1:n-1
        for j=i+1:n
            G[i,j]=G[j,i]= 1 / (2*pi) * (1 + A * cos(u[i] - u[j]))
        end
    end

    return G
end

function palacios(du,u,p,t)
    global omega, G, alpha, n
    for i = 1:n
        du[i] = omega[i] - ((sum(j->(G[i,j] * sin(u[i]-u[j]-alpha)),1:n))) 
    end
    
end


open("heat_phases.txt", "w") do io
    write(io, "A\talpha\ty\n")
end

global n = 256                                 
global omega = 0.01*ones(n)   
global B = 0

while B <= 1
    global B
    global beta = 0
    while beta <= 0.3
        
        global beta
        
        rng = MersenneTwister(58959);
        theta0=[6*(rand(rng)-0.5)*exp(-30*(x/float(n)-0.5)*(x/float(n)-0.5))  for x in range(1,stop=n)]
        x0 = Array(range(0,stop=2*pi,length=n))
        A = B
        global G = regular(x0, n, A);
        global alpha = pi/2 - beta

        niter=500               
        dt = 0.025                  
        ti = 0.0; 
        tf = niter*dt
        tspan = (ti, tf)             

        #system integration
        prob = ODEProblem(palacios,theta0,tspan)

        sol = solve(prob, Tsit5(), saveat=dt) ; #reltol=1e-6, 
        
        beta += 0.01
        

        s =join(sol[:, niter], ';')
        open("heat_phases.txt", "a") do io
            write(io, "$(A)\t$(alpha)\t$(s)\n")
        end
    
    end
    println(B)
    B += 0.01
end
           
        # Get if the persenatge of mode count?
        
        # Print on a file A, alpha, density