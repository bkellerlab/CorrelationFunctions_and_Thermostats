using Chemfiles
using LinearAlgebra
using PyCall
using JLD2
using LaTeXStrings


flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
GC.gc()
plt=pyimport("matplotlib.pyplot")
gs=pyimport("matplotlib.gridspec")
pygui()


rNVE=[1,6.25,0.0876]
rNH=[2,6.48,0.122]
rBe=[3,6.52,0.177]

rvr_001=[4.5,6.44,0.0601]
rvr_01=[5.5,6.63,0.187]
rvr_1=[6.5,6.48,0.109]
rvr_10=[7.5,6.41,0.105]
rvr_100=[8.5,6.38,0.23]
rsd_001=[10,1.36,0.00971]
rsd_01=[11,3.37,0.0169]
rsd_1=[12,5.7,0.0541]
rsd_10=[13,6.35,0.131]
rsd_100=[14,6.41,0.216]

xtb_NVE=[18,4.41,0.173]

xtb_vr_0001=[19.5,4.44,0.165]
xtb_vr_001=[20.5,4.17,0.164]
xtb_vr_01=[21.5,4.39,0.162]
xtb_vr_1=[22.5,4.69,0.141]

xtb_sd_0001=[24,0.09,0.0024]
xtb_sd_001=[25,0.6,0.018]
xtb_sd_01=[26,2.71,0.123]
xtb_sd_1=[27,3.98,0.123]

be=load("VACF_TIP3P_JLD2/be_vacf.jld2")
nh=load("VACF_TIP3P_JLD2/nh_vacf.jld2")
nve=load("VACF_TIP3P_JLD2/tip3p_vacf.jld2")
vr_001=load("VACF_TIP3P_JLD2/vr_001_vacf.jld2")
vr_01=load("VACF_TIP3P_JLD2/vr_01_vacf.jld2")
vr_1=load("VACF_TIP3P_JLD2/vr_1_vacf.jld2")
vr_10=load("VACF_TIP3P_JLD2/vr_10_vacf.jld2")
vr_100=load("VACF_TIP3P_JLD2/vr_100_vacf.jld2")

sd_001=load("VACF_TIP3P_JLD2/sd_001_vacf.jld2")
sd_01=load("VACF_TIP3P_JLD2/sd_01_vacf.jld2")
sd_1=load("VACF_TIP3P_JLD2/sd_1_vacf.jld2")
sd_10=load("VACF_TIP3P_JLD2/sd_10_vacf.jld2")
sd_100=load("VACF_TIP3P_JLD2/sd_100_vacf.jld2")

times=sd_10["times"]*1E+12

A=Dict("width_ratios"=>[6],"height_ratios"=>[2.5,4])

fig, ax = plt.subplots(2,1,figsize=(6, 7),gridspec_kw=A)
ax[1,1].plot(times,vr_001["VACF_N"],label=L"VR $τ_T$=0.01 ps",color="#D32F2F",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_01["VACF_N"],label=L"VR $τ_T$=0.1 ps",color="#800080",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_1["VACF_N"],label=L"VR $τ_T$=1 ps",color="#00509E",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_10["VACF_N"],label=L"VR $τ_T$=10 ps",color="#006400",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_100["VACF_N"],label=L"VR $τ_T$=100ps",color="#FFCC00",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,sd_001["VACF_N"],label=L"GSD $τ_T$=0.01 ps",color="#D32F2F",linewidth=1)
ax[1,1].plot(times,sd_01["VACF_N"],label=L"GSD $τ_T$=0.1 ps",color="#800080",linewidth=1)
ax[1,1].plot(times,sd_1["VACF_N"],label=L"GSD $τ_T$=1 ps",color="#00509E",linewidth=1)
ax[1,1].plot(times,sd_10["VACF_N"],label=L"GSD $τ_T$=10 ps",color="#006400",linewidth=1)
ax[1,1].plot(times,sd_100["VACF_N"],label=L"GSD $τ_T$=100ps",color="#FFCC00",linewidth=1)
ax[1,1].set_xlabel("time in ps", fontsize=14)
ax[1,1].set_ylabel("VACF",fontsize=14)
ax[1,1].tick_params(axis="both",which="major",labelsize=14)
ax[1,1].set_xlim(0,1)
ax[1,1].set_ylim(-0.05,1)
ax[1,1].set_yticks([0,0.5,1])
ax[1,1].set_xticks([0,0.2,0.4,0.6,0.8,1])
ax[1,1].set_xticklabels(["0","0.2","0.4","0.6","0.8","1"])
ax[1,1].legend(fontsize=14,ncol=2,borderpad=0.3,handletextpad=0.3,handlelength=1)
ax[1,1].set_title("a)",fontsize=18,y=1,pad=-30,x=-0.11)

ax[2,1].bar(rNVE[1],rNVE[2],yerr=rNVE[3],label="NVE",color="#000000",capsize=2,width=1)
ax[2,1].bar(rNH[1],rNH[2],yerr=rNH[3],label=L"NH $τ_T^{NH}$= 4ps",color="#003366",capsize=2,width=1)
ax[2,1].bar(rBe[1],rBe[2],yerr=rBe[3],label=L"BE $τ_T$= 1ps",color="#008B8B",width=1,capsize=2)

ax[2,1].bar(xtb_vr_0001[1],xtb_vr_0001[2],yerr=xtb_vr_0001[3],label=L"VR $τ_T$= 0.001ps",color="#FF8C00",width=1,capsize=2, hatch="//",alpha=0.99)
ax[2,1].bar(rvr_001[1],rvr_001[2],yerr=rvr_001[3],label=L"VR $τ_T$= 0.01ps",color="#D32F2F",width=1,capsize=2, hatch="//", alpha=0.99)
ax[2,1].bar(rvr_01[1],rvr_01[2],yerr=rvr_01[3],label=L"VR $τ_T$= 0.1ps",color="#800080",width=1,capsize=2, hatch="//", alpha=0.99)
ax[2,1].bar(rvr_1[1],rvr_1[2],yerr=rvr_1[3],label=L"VR $τ_T$= 1ps",color="#00509E",width=1,capsize=2, hatch="//", alpha=0.99)
ax[2,1].bar(rvr_10[1],rvr_10[2],yerr=rvr_10[3],label=L"VR $τ_T$= 10ps",color="#006400",width=1,capsize=2, hatch="//", alpha=0.99)
ax[2,1].bar(rvr_100[1],rvr_100[2],yerr=rvr_100[3],label=L"VR $τ_T$= 100ps",color="#FFCC00",width=1,capsize=2, hatch="//", alpha=0.99)

ax[2,1].bar(xtb_sd_0001[1],xtb_sd_0001[2],yerr=xtb_sd_0001[3],label=L"SD $τ_T$= 0.001ps",color="#FF8C00",width=1,capsize=2)
ax[2,1].bar(rsd_001[1],rsd_001[2],yerr=rsd_001[3],label=L"SD $τ_T$= 0.01ps",color="#D32F2F",width=1,capsize=2)
ax[2,1].bar(rsd_01[1],rsd_01[2],yerr=rsd_01[3],label=L"SD $τ_T$= 0.1ps",color="#800080",width=1,capsize=2)
ax[2,1].bar(rsd_1[1],rsd_1[2],yerr=rsd_1[3],label=L"SD $τ_T$= 1ps",color="#00509E",width=1,capsize=2)
ax[2,1].bar(rsd_10[1],rsd_10[2],yerr=rsd_10[3],label=L"SD $τ_T$= 10ps",color="#006400",width=1,capsize=2)
ax[2,1].bar(rsd_100[1],rsd_100[2],yerr=rsd_100[3],label=L"SD $τ_T$= 100ps",color="#FFCC00",width=1,capsize=2)

ax[2,1].bar(xtb_NVE[1],xtb_NVE[2],yerr=xtb_vr_001[3],color="#000000",width=1,capsize=2)
ax[2,1].bar(xtb_vr_001[1],xtb_vr_001[2],yerr=xtb_vr_001[3],color="#D32F2F",width=1,capsize=2, hatch="//", alpha=0.99)
ax[2,1].bar(xtb_vr_01[1],xtb_vr_01[2],yerr=xtb_vr_01[3],color="#800080",width=1,capsize=2, hatch="//", alpha=0.99)
ax[2,1].bar(xtb_vr_1[1],xtb_vr_1[2],yerr=xtb_vr_1[3],color="#00509E",width=1,capsize=2, hatch="//", alpha=0.99)
ax[2,1].bar(xtb_sd_001[1],xtb_sd_001[2],yerr=xtb_sd_001[3],color="#D32F2F",width=1,capsize=2)
ax[2,1].bar(xtb_sd_01[1],xtb_sd_01[2],yerr=xtb_sd_01[3],color="#800080",width=1,capsize=2)
ax[2,1].bar(xtb_sd_1[1],xtb_sd_1[2],yerr=xtb_sd_1[3],color="#00509E",width=1,capsize=2)

ax[2,1].set_ylabel(L"D in $10^{-5}$ cm/s",fontsize=14)
ax[2,1].tick_params(axis="both",which="major",labelsize=14)
ax[2,1].set_xticks([7,22])
ax[2,1].set_xticklabels(["TIP3P","xTB"])
ax[2,1].set_ylim([0,10])
ax[2,1].set_xlim([0,28])
ax[2,1].legend(fontsize=11,loc="upper left",labelspacing = 0.1,ncol=3, borderpad=0.3,handletextpad=0.3,handlelength=1.5,columnspacing=0.7)
ax[2,1].set_title("b)",fontsize=18,y=1,pad=-30,x=-0.11)
plt.tight_layout()
plt.savefig("Diffusion_combined.pdf")
plt.close()
