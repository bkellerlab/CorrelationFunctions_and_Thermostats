using DelimitedFiles
using PyCall
using Statistics
using Latexify
using StatsBase
using LinearAlgebra
using Distributions
using LaTeXStrings

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
GC.gc()
plt = pyimport("matplotlib.pyplot")
pygui()

sd_002=readdlm("sd_002.xvg",skipstart=25)
sd_02=readdlm("sd_02.xvg",skipstart=25)
sd_2=readdlm("sd_2.xvg",skipstart=25)
sd_20=readdlm("sd_20.xvg",skipstart=25)
sd_100=readdlm("sd_100.xvg",skipstart=25)
vr_001=readdlm("vr_001.xvg",skipstart=25)
vr_01=readdlm("vr_01.xvg",skipstart=25)
vr_1=readdlm("vr_1.xvg",skipstart=25)
vr_10=readdlm("vr_10.xvg",skipstart=25)
vr_100=readdlm("vr_100.xvg",skipstart=25)
be=readdlm("be_1.xvg",skipstart=25)
be_vv=readdlm("be_1_md-vv.xvg",skipstart=25)
nh=readdlm("nh_2.xvg",skipstart=25)
nve=readdlm("nve2.xvg",skipstart=25)


range=collect(-33.7:0.05:-31.6)
range2=collect(-33.0:0.01:-32.5)
range3=collect(-32.8:0.001:-32.5)
vr_1_h=normalize(fit(Histogram,vr_1[10001:end,3]/884,range),mode=:pdf)
be_h=normalize(fit(Histogram,be[10001:end,3]/884,range2),mode=:pdf)
nh_h=normalize(fit(Histogram,nh[10001:end,3]/884,range),mode=:pdf)
sd_2_h=normalize(fit(Histogram,sd_2[10001:end,3]/884,range),mode=:pdf)
sd_100_h=normalize(fit(Histogram,sd_100[10001:end,3]/884,range),mode=:pdf)
nve_h=normalize(fit(Histogram,nve[10001:end,3]/884,range3),mode=:pdf)


b=fit(Normal,nh[10001:end,3]/884)


T=300                   #K
Cv=4.14                 #kJ/kgK
R=8.314/1000            #kJ/molK
M=0.018                 #kg/mol

X=-32.684               #kJ/mol
s=sqrt(R*Cv*T^2*M/884)           #kJ/molK*kJ/kgK*K^2*kg/mol=(kJ/mol)^2
function gauss(x)       
    return 1/sqrt(2*3.14159*s^2)*exp(-(x-X)^2/(2*s^2))
end

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(range[1:end-1],vr_1_h.weights./100,label="Velocity-Rescale" ,color="#00509E",linestyle="dashed",linewidth=1.5)
ax.plot(range2[1:end-1],be_h.weights./100,label="Berendsen",color="#008B8B",linewidth=2.0)
ax.plot(range[1:end-1],nh_h.weights./100,label="Nose-Hoover",color="#003366",linewidth=1.0)
ax.plot(range[1:end-1],sd_2_h.weights./100,label=L"GSD $\tau_T$=2 ps",color="#00509E",linewidth=1.0)
ax.plot(range[1:end-1],sd_100_h.weights./100,label=L"GSD $\tau_T$=100 ps",color="#FFCC00",linewidth=1.0)
ax.plot(range3[1:end-1].+0.12,nve_h.weights./100,label="NVE",color="#000000",linewidth=1.0)
ax.plot(range,gauss.(range)./100,label="NVT",color="red",linestyle="-.",linewidth=1.5)
ax.set_xlabel("Total Energy in kJ/mol", fontsize=14)
ax.set_ylabel("P",fontsize=14)
ax.legend(loc="upper right",fontsize=14, borderpad=0.3,handletextpad=0.3,handlelength=1.5,framealpha=1)
ax.set_ylim(0,0.095)
ax.tick_params(axis="both",which="major",labelsize=14)
plt.tight_layout()
plt.savefig("EnergyDistribution.pdf")
plt.close()

#########################################################################################################
range=collect(-33.7:0.05:-31.6)

vr_100_h=normalize(fit(Histogram,vr_100[:,3]./884,range),mode=:pdf)
vr_10_h=normalize(fit(Histogram,vr_10[:,3]./884,range),mode=:pdf)
vr_1_h=normalize(fit(Histogram,vr_1[:,3]./884,range),mode=:pdf)
vr_01_h=normalize(fit(Histogram,vr_01[:,3]./884,range),mode=:pdf)
vr_001_h=normalize(fit(Histogram,vr_001[:,3]./884,range),mode=:pdf)

sd_100_h=normalize(fit(Histogram,sd_100[:,3]./884,range),mode=:pdf)
sd_20_h=normalize(fit(Histogram,sd_20[:,3]./884,range),mode=:pdf)
sd_2_h=normalize(fit(Histogram,sd_2[:,3]./884,range),mode=:pdf)
sd_02_h=normalize(fit(Histogram,sd_02[:,3]./884,range),mode=:pdf)
sd_002_h=normalize(fit(Histogram,sd_002[:,3]./884,range),mode=:pdf)

fig, ax = plt.subplots()
ax.plot(range[1:end-1],vr_001_h.weights./100,label=L"VR $\tau_T$=0.01 ps",color="#D32F2F",linestyle="dashed",linewidth=1.0)
ax.plot(range[1:end-1],vr_01_h.weights./100,label=L"VR $\tau_T$=0.1 ps",color="#800080",linestyle="dashed",linewidth=1.0)
ax.plot(range[1:end-1],vr_1_h.weights./100,label=L"VR $\tau_T$=1.0 ps",color="#00509E",linestyle="dashed",linewidth=1.0)
ax.plot(range[1:end-1],vr_10_h.weights./100,label=L"VR $\tau_T$=10 ps",color="#006400",linestyle="dashed",linewidth=1.0)
ax.plot(range[1:end-1],vr_100_h.weights./100,label=L"VR $\tau_T$=100 ps",color="#FFCC00",linestyle="dashed",linewidth=1.0)
ax.set_xlabel("Total Energy in kJ/mol", fontsize=14)
ax.set_ylabel("P",fontsize=14)
ax.legend(loc="upper right",fontsize=14, borderpad=0.3,handletextpad=0.6,handlelength=1.2,framealpha=0.95)
ax.set_ylim(0,0.020)
ax.tick_params(axis="both",which="major",labelsize=14)
plt.tight_layout()
#plt.show()
plt.savefig("EnergyDistribution_VR.pdf")
plt.close()

fig, ax = plt.subplots()
ax.plot(range[1:end-1],sd_002_h.weights./100,label=L"VR $\tau_T$=0.02 ps",color="#D32F2F",linewidth=1.0)
ax.plot(range[1:end-1],sd_02_h.weights./100,label=L"VR $\tau_T$=0.2 ps",color="#800080",linewidth=1.0)
ax.plot(range[1:end-1],sd_2_h.weights./100,label=L"VR $\tau_T$=2.0 ps",color="#00509E",linewidth=1.0)
ax.plot(range[1:end-1],sd_20_h.weights./100,label=L"VR $\tau_T$=20 ps",color="#006400",linewidth=1.0)
ax.plot(range[1:end-1],sd_100_h.weights./100,label=L"VR $\tau_T$=100 ps",color="#FFCC00",linewidth=1.0)
ax.set_xlabel("Total Energy in kJ/mol", fontsize=14)
ax.set_ylabel("P",fontsize=14)
ax.legend(loc="upper right",fontsize=14, borderpad=0.3,handletextpad=0.6,handlelength=1.2,framealpha=1)
ax.set_ylim(0,0.020)
ax.tick_params(axis="both",which="major",labelsize=14)
plt.tight_layout()
#plt.show()
plt.savefig("EnergyDistribution_SD.pdf")
plt.close()

#################################


range4=collect(-15182:0.1:-15176)
sd_1=[readdlm("Frederick/E_tot-SD-tau=1ps.dat")...].*2625.5./128
sd_01=[readdlm("Frederick/E_tot-SD-tau=0.1ps.dat")...].*2625.5./128
sd_001=[readdlm("Frederick/E_tot-SD-tau=0.01ps.dat")...].*2625.5./128
sd_0001=[readdlm("Frederick/E_tot-SD-tau=0.001ps.dat")...].*2625.5./128


sd_1_h=normalize(fit(Histogram,sd_1,range4),mode=:pdf)
sd_01_h=normalize(fit(Histogram,sd_01,range4),mode=:pdf)
sd_001_h=normalize(fit(Histogram,sd_001,range4),mode=:pdf)
sd_0001_h=normalize(fit(Histogram,sd_0001,range4),mode=:pdf)

fig, ax = plt.subplots()
ax.plot(range4[1:end-1],sd_0001_h.weights./100,label=L"SD $\tau_T$=0.001 ps",color="#FF8C00")
ax.plot(range4[1:end-1],sd_001_h.weights./100,label=L"SD $\tau_T$=0.01 ps",color="#D32F2F")
ax.plot(range4[1:end-1],sd_01_h.weights./100,label=L"SD $\tau_T$=0.1 ps",color="#800080")
ax.plot(range4[1:end-1],sd_1_h.weights./100,label=L"SD $\tau_T$=1.0 ps",color="#00509E")
ax.set_xlabel("Total Energy in kJ/mol", fontsize=16)
ax.set_ylabel("P",fontsize=16)
ax.legend(loc="upper left",fontsize=16)
ax.set_ylim(0,0.02)
ax.tick_params(axis="both",which="major",labelsize=14)
plt.tight_layout()
plt.savefig("EnergyDistribution_xTB_SD.pdf")
plt.close()


vr_1=[readdlm("Frederick/E_tot-VR-tau=1ps.dat")...].*2625.5./128
vr_01=[readdlm("Frederick/E_tot-VR-tau=0.1ps.dat")...].*2625.5./128
vr_001=[readdlm("Frederick/E_tot-VR-tau=0.01ps.dat")...].*2625.5./128
vr_0001=[readdlm("Frederick/E_tot-VR-tau=0.001ps.dat")...].*2625.5./128

vr_1_h=normalize(fit(Histogram,vr_1,range4),mode=:pdf)
vr_01_h=normalize(fit(Histogram,vr_01,range4),mode=:pdf)
vr_001_h=normalize(fit(Histogram,vr_001,range4),mode=:pdf)
vr_0001_h=normalize(fit(Histogram,vr_0001,range4),mode=:pdf)

fig, ax = plt.subplots()
ax.plot(range4[1:end-1],vr_0001_h.weights./100,label=L"VR $\tau_T$=0.001 ps",color="#FF8C00",linestyle="dashed")
ax.plot(range4[1:end-1],vr_001_h.weights./100,label=L"VR $\tau_T$=0.01 ps",color="#D32F2F",linestyle="dashed")
ax.plot(range4[1:end-1],vr_01_h.weights./100,label=L"VR $\tau_T$=0.1 ps",color="#800080",linestyle="dashed")
ax.plot(range4[1:end-1],vr_1_h.weights./100,label=L"VR $\tau_T$=1.0 ps",color="#00509E",linestyle="dashed")
ax.set_xlabel("Total Energy in kJ/mol", fontsize=16)
ax.set_ylabel("P",fontsize=16)
ax.legend(loc="upper right",fontsize=16)
ax.set_ylim(0,0.02)
ax.tick_params(axis="both",which="major",labelsize=13)
plt.tight_layout()
plt.savefig("EnergyDistribution_xTB_VR.pdf")
plt.close()


println("END")