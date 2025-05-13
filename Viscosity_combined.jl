using Chemfiles
using LinearAlgebra
using PyCall
using JLD2
using LaTeXStrings

using StatsBase
using NumericalIntegration

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
GC.gc()
plt=pyimport("matplotlib.pyplot")
pygui()

s=["Pentane","Water","Aniline","Glycerine"]
EXP=[0.24,1,4.40,1426]
NH_mean=[0.107,0.339,1.92,1710]
NH_STD=[0.0266,0.0659,0.849,1680]
BE_mean=[0.105,0.319,1.83,1880]
BE_STD=[0.0266,0.0389,0.636,1440]
VR_1_mean=[0.116,0.334,2.0,1600]
VR_1_std=[0.0256,0.0468,0.886,1450]
SD_001_mean=[1.69,2.62,23.4,2220]
SD_001_std=[0.665,0.478,17.4,1410]
SD_01_mean=[0.229,0.583,6.1,1170]
SD_01_std=[0.0701,0.109,3.78,905]
SD_1_mean=[0.139,0.369,2.2,992]
SD_1_std=[0.0476,0.0633,1.42,713]
SD_10_mean=[0.114,0.355,1.75,1910]
SD_10_std=[0.0362,0.0667,0.761,999]
SD_100_mean=[0.101,0.358,2.74,2180]
SD_100_std=[0.0265,0.0535,1.71,1780]

v=[10,21,32,43]


vr_001=load("PACF_JLD2/TIP3P_JLD2/vr_001_pacf.jld2")
vr_01=load("PACF_JLD2/TIP3P_JLD2/vr_01_pacf.jld2")
vr_1=load("PACF_JLD2/TIP3P_JLD2/vr_1_pacf.jld2")
vr_10=load("PACF_JLD2/TIP3P_JLD2/vr_10_pacf.jld2")
vr_100=load("PACF_JLD2/TIP3P_JLD2/vr_100_pacf.jld2")

sd_001=load("PACF_JLD2/TIP3P_JLD2/sd_001_pacf.jld2")
sd_01=load("PACF_JLD2/TIP3P_JLD2/sd_01_pacf.jld2")
sd_1=load("PACF_JLD2/TIP3P_JLD2/sd_1_pacf.jld2")
sd_10=load("PACF_JLD2/TIP3P_JLD2/sd_10_pacf.jld2")
sd_100=load("PACF_JLD2/TIP3P_JLD2/sd_100_pacf.jld2")

times=sd_10["times"]*1E+12
times2=sd_10["times2"]*1E+12


name2="Anilin"
pvr_1=load("PACF_JLD2/"*name2*"_JLD2/vr_1_pacf.jld2")
psd_001=load("PACF_JLD2/"*name2*"_JLD2/sd_001_pacf.jld2")
psd_01=load("PACF_JLD2/"*name2*"_JLD2/sd_01_pacf.jld2")
psd_1=load("PACF_JLD2/"*name2*"_JLD2/sd_1_pacf.jld2")
psd_10=load("PACF_JLD2/"*name2*"_JLD2/sd_10_pacf.jld2")
psd_100=load("PACF_JLD2/"*name2*"_JLD2/sd_100_pacf.jld2")
ptimes=sd_10["times"]*1E+12
ptimes2=sd_10["times2"]*1E+12

A=Dict("width_ratios"=>[6],"height_ratios"=>[3,3,1,4])

fig, ax = plt.subplots(4,1,figsize=(6, 11),gridspec_kw=A)

ax[1,1].plot(times,vr_001["PACF_N"],label=L"VR $τ_T$=0.01 ps",color="#D32F2F",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_01["PACF_N"],label=L"VR $τ_T$=0.1 ps",color="#800080",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_1["PACF_N"],label=L"VR $τ_T$=1 ps",color="#00509E",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_10["PACF_N"],label=L"VR $τ_T$=10 ps",color="#006400",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,vr_100["PACF_N"],label=L"VR $τ_T$=100ps",color="#FFCC00",linestyle="dashed",linewidth=2)
ax[1,1].plot(times,sd_001["PACF_N"],label=L"GSD $τ_T$=0.01 ps",color="#D32F2F",linewidth=1)
ax[1,1].plot(times,sd_01["PACF_N"],label=L"GSD $τ_T$=0.1 ps",color="#800080",linewidth=1)
ax[1,1].plot(times,sd_1["PACF_N"],label=L"GSD $τ_T$=1 ps",color="#00509E",linewidth=1)
ax[1,1].plot(times,sd_10["PACF_N"],label=L"GSD $τ_T$=10 ps",color="#006400",linewidth=1)
ax[1,1].plot(times,sd_100["PACF_N"],label=L"GSD $τ_T$=100ps",color="#FFCC00",linewidth=1)
ax[1,1].set_xlabel("time in ps", fontsize=14)
ax[1,1].set_ylabel("PACF",fontsize=14)
ax[1,1].tick_params(axis="both",which="major",labelsize=14)
ax[1,1].set_xlim(0,1)
ax[1,1].set_ylim(-0.1,1)
ax[1,1].set_xticks([0,0.2,0.4,0.6,0.8,1])
ax[1,1].set_xticklabels(["0","0.2","0.4","0.6","0.8","1"])
ax[1,1].legend(loc="upper right",fontsize=13,ncol=2,borderpad=0.3,handletextpad=0.3,handlelength=1)
ax[1,1].set_title("a)",fontsize=18,y=1,pad=-25,x=-0.14)

ax[2,1].plot(ptimes,pvr_1["PACF_N"],label=L"VR $τ_T$=1 ps",color="#00509E",linestyle="dashed",linewidth=2)
ax[2,1].plot(ptimes,psd_001["PACF_N"],label=L"GSD $τ_T$=0.01 ps",color="#D32F2F",linewidth=1)
ax[2,1].plot(ptimes,psd_01["PACF_N"],label=L"GSD $τ_T$=0.1 ps",color="#800080",linewidth=1)
ax[2,1].plot(ptimes,psd_1["PACF_N"],label=L"GSD $τ_T$=1 ps",color="#00509E",linewidth=1)
ax[2,1].plot(ptimes,psd_10["PACF_N"],label=L"GSD $τ_T$=10 ps",color="#006400",linewidth=1)
ax[2,1].plot(ptimes[1:1000],psd_100["PACF_N"][1:1000],label=L"GSD $τ_T$=100ps",color="#FFCC00",linewidth=1)
ax[2,1].set_xlabel("time in ps", fontsize=14)
ax[2,1].set_ylabel("PACF",fontsize=14)
ax[2,1].tick_params(axis="both",which="major",labelsize=14)
ax[2,1].set_xticks([0,0.2,0.4,0.6,0.8,1])
ax[2,1].set_xticklabels(["0","0.2","0.4","0.6","0.8","1"])
ax[2,1].set_xlim(0,1)
ax[2,1].set_ylim(-0.4,1)
ax[2,1].legend(loc="upper right",fontsize=14,ncol=2,borderpad=0.3,handletextpad=0.3,handlelength=1)
ax[2,1].set_title("b)",fontsize=18,y=1,pad=-30,x=-0.13)

ax[3,1].plot(ptimes,psd_001["PACF_N"],label=L"GSD $τ_T$=0.01 ps",color="#D32F2F",linewidth=1)
ax[3,1].plot(ptimes,psd_1["PACF_N"],label=L"GSD $τ_T$=1 ps",color="#00509E",linewidth=1)
ax[3,1].set_xlabel("time in ps", fontsize=14)
ax[3,1].set_ylabel("PACF",fontsize=14)
ax[3,1].tick_params(axis="both",which="major",labelsize=14)
ax[3,1].set_xlim(0,1)
ax[3,1].set_xticks([0,25,50,75,100,125,150,175])
ax[3,1].set_xticklabels(["0","25","50","75","100","125","150","175"])
ax[3,1].legend(loc="upper right",fontsize=14,ncol=2,borderpad=0.3,handletextpad=0.3,handlelength=1)
ax[3,1].set_title("c)",fontsize=18,y=1,pad=0,x=-0.14)
ax[3,1].set_xlim(0,200)
ax[3,1].set_ylim(0,0.02)

ax[4,1].bar(v.-5,EXP,label="Experimental",color="grey",width=1)
ax[4,1].bar(v.-4,NH_mean,yerr=NH_STD,label="NH τ_c= 4ps",color="#003366",capsize=2,width=1)
ax[4,1].bar(v.-3,BE_mean,yerr=BE_STD,label="BE τ= 1ps",color="#008B8B",width=1,capsize=2)
ax[4,1].bar(v.-2,VR_1_mean,yerr=VR_1_std,label="VR τ= 1ps",color="#00509E",width=1,capsize=2, hatch="//",alpha=0.99)
ax[4,1].bar(v.-0.5,SD_001_mean,yerr=SD_001_std,label="SD τ= 0.01ps",color="#D32F2F",width=1,capsize=2)
ax[4,1].bar(v.+0.5,SD_01_mean,yerr=SD_01_std,label="SD τ= 0.1ps",color="#800080",width=1,capsize=2)
ax[4,1].bar(v.+1.5,SD_1_mean,yerr=SD_1_std,label="SD τ= 1ps",color="#00509E",width=1,capsize=2)
ax[4,1].bar(v.+2.5,SD_10_mean,yerr=SD_10_std,label="SD τ= 10ps",color="#006400",width=1,capsize=2)
ax[4,1].bar(v.+3.5,SD_100_mean,yerr=SD_100_std,label="SD τ= 100ps",color="#FFCC00",width=1,capsize=2)
ax[4,1].set_ylabel("viscosity η in mPas",fontsize=14)
ax[4,1].tick_params(axis="both",which="major",labelsize=14)
ax[4,1].set_yscale("log")
ax[4,1].set_xticks([9,21,31,43])
ax[4,1].set_xticklabels(["Pentane","Water","Aniline","Glycerol"])
ax[4,1].set_ylim([0.06,10000])
ax[4,1].set_xlim([4,47.5])
ax[4,1].set_yticks([0.1,1,10,100,1000])
ax[4,1].set_yticklabels(["0.1","1","10","100","1000"])
ax[4,1].legend(fontsize=13,loc="upper left",labelspacing = 0.1,ncol=2, borderpad=0.3,handletextpad=0.3,handlelength=1.5,columnspacing=0.7)
ax[4,1].set_title("d)",fontsize=18,y=1,pad=-10,x=-0.14)

plt.tight_layout()
plt.savefig("viscosity_result.pdf")
plt.close()

println("End")
