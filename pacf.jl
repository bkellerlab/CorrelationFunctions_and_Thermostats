using Chemfiles
using LinearAlgebra
using PyCall
using JLD2
using FFTW
using StatsBase
using DelimitedFiles
using NumericalIntegration

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
GC.gc()
plt=pyimport("matplotlib.pyplot")
pygui()


V=27E-27
k=1.38E-23
T=300
XX=V/(k*T)

const start_water=1804
const m_O=15.999
const m_H=1.00784
const m_H20=m_O+2m_H

#calculates the center of mass positionn and center of mass velocietes for Water molecules
function center_of_mass(values)
    com=(m_O*values[:,1:3:end] .+m_H*values[:,2:3:end].+m_H*values[:,2:3:end])/m_H20
    return com
end

#Transforms two points into an vector that conforms to periodic boundary conditions
function pbc_dir_vector(A,B,box,r_box)
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    return dx
end

#Mean position between two points that conforms to periodic boundary conditions
function pbc_mean(A,B,box,r_box)
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    A=A-0.5*dx
    return A
end

#Calculates the velocity autocorraletion function
function ACF(Values)
    N=length(Values)
    X=zeros(2*N)
    X[1:N]=Values
    X=fft(X)
    X=X.*conj.(X)
    X=ifft(X)
    X=real.(X[1:N])
    return X./N
end


function PACF(name,header::Int64,batches::Int64)
    data=readdlm("$name.xvg",skipstart=header,Float64)
    times=data[:,1]*1E-12
    data[:,2:end]=data[:,2:end]*10^5
    t1=floor(Int64,length(times)/batches)
    PACF=ACF(data[:,2])./3 .+ACF(data[:,3])./3 .+ACF(data[:,4])./3
    pacf=zeros(t1,batches)
    for i in 1:batches
        pacf[:,i]=ACF(data[(1+(i-1)t1):(i*t1),2])./3 .+ACF(data[(1+(i-1)t1):(i*t1),4])./3 .+ACF(data[(1+(i-1)t1):(i*t1),4])./3 
    end
    pacf_mean=[mean(pacf,dims=2)...]
    pacf_mean=pacf_mean./pacf_mean[1]
    pacf_std=[std(pacf,dims=2)...]

    
    return PACF,pacf,pacf_mean,pacf_std,times
end
@time nh,nh2,nh_mean,nh_std,times=PACF("nh_pacf2",26,100)
@time vr,vr2,vr_mean,vr_std,times=PACF("vr_pacf2",26,100)
@time be,be2,be_mean,be_std,times=PACF("be_pacf2",26,100)
@time nve,nve2,nve_mean,nve_std,times=PACF("nve_pacf2",26,100)

fig, ax = plt.subplots()
ax.plot(times[1:201],vr_mean[1:201],label="Velocity-Rescale",color="red",linestyle="--")
ax.plot(times[1:201],be_mean[1:201],label="Berendsen",color="green",linestyle="--")
ax.plot(times[1:201],nh_mean[1:201],label="Nose-Hoover",color="blue",linestyle="--")
ax.plot(times[1:201],nve_mean[1:201],label="NVE",color="black",linestyle="--")
ax.set_xlabel("time in ps", fontsize=16)
ax.set_ylabel("PACF",fontsize=16)
ax.tick_params(axis="both",which="major",labelsize=16)
plt.tight_layout()
plt.legend(fontsize=16)
plt.savefig("PACF3.pdf")
plt.close()
#plt.show()

V=27E-27
k=1.381E-23
T=300
XX=V/(k*T)
t=times[1:40000]
t2=times[1:40000]*1E+12

nh_c_mean=[mean(cumul_integrate(t,nh2,dims=2),dims=2)...]*1000*XX
be_c_mean=[mean(cumul_integrate(t,be2,dims=2),dims=2)...]*1000*XX
vr_c_mean=[mean(cumul_integrate(t,vr2,dims=2),dims=2)...]*1000*XX
nve_c_mean=[mean(cumul_integrate(t,nve2,dims=2),dims=2)...]*1000*XX

nh_c_std=[std(cumul_integrate(t,nh2,dims=2),dims=2)...]*1000*XX
be_c_std=[std(cumul_integrate(t,be2,dims=2),dims=2)...]*1000*XX
vr_c_std=[std(cumul_integrate(t,vr2,dims=2),dims=2)...]*1000*XX
nve_c_std=[std(cumul_integrate(t,nve2,dims=2),dims=2)...]*1000*XX


fig, ax = plt.subplots()
ax.plot(t2,vr_c_mean,label="Velocity-Rescale",color="red")
ax.fill_between(t2,vr_c_mean.-vr_c_std,vr_c_mean.+vr_c_std,color="red",alpha=0.2)
ax.plot(t2,be_c_mean,label="Berendsen",color="green")
ax.fill_between(t2,be_c_mean.-be_c_std,be_c_mean.+be_c_std,color="green",alpha=0.2)
ax.plot(t2,nh_c_mean,label="Nose-Hoover",color="blue")
ax.fill_between(t2,nh_c_mean.-nh_c_std,nh_c_mean.+nh_c_std,color="blue",alpha=0.2)
ax.plot(t2,nve_c_mean,label="NVE",color="black")
ax.fill_between(t2,nve_c_mean.-nve_c_std,nve_c_mean.+nve_c_std,color="black",alpha=0.2)
ax.set_xlabel("time in ps", fontsize=16)
ax.set_ylabel("viscosity in mPas",fontsize=16)
ax.tick_params(axis="both",which="major",labelsize=16)
ax.set_xlim(0,20)
plt.tight_layout()
plt.legend(fontsize=16)
plt.savefig("PACF4.pdf")
plt.show()
