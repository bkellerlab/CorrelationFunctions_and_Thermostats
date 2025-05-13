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

T=1000
time_v=load("VACF_TIP3P_JLD2/vr_1_vacf.jld2")["times2"][1:T]
nve_v=load("VACF_TIP3P_JLD2/vr_1_vacf.jld2")["vacf_mean"][1:T]
sd_001_v=load("VACF_TIP3P_JLD2/sd_001_vacf.jld2")["vacf_mean"][1:T]
sd_01_v=load("VACF_TIP3P_JLD2/sd_01_vacf.jld2")["vacf_mean"][1:T]
sd_1_v=load("VACF_TIP3P_JLD2/sd_1_vacf.jld2")["vacf_mean"][1:T]
sd_10_v=load("VACF_TIP3P_JLD2/sd_10_vacf.jld2")["vacf_mean"][1:T]
sd_100_v=load("VACF_TIP3P_JLD2/sd_100_vacf.jld2")["vacf_mean"][1:T]
nve_v_int=integrate(time_v,abs.(nve_v))

v_001=integrate(time_v,abs.(sd_001_v.-nve_v))/nve_v_int
v_01=integrate(time_v,abs.(sd_01_v.-nve_v))/nve_v_int
v_1=integrate(time_v,abs.(sd_1_v.-nve_v))/nve_v_int
v_10=integrate(time_v,abs.(sd_10_v.-nve_v))/nve_v_int
v_100=integrate(time_v,abs.(sd_100_v.-nve_v))/nve_v_int
tip3p_sd_v=[v_001,v_01,v_1,v_10,v_100]


T=200
time_p=load("PACF_JLD2/TIP3P_JLD2/vr_1_pacf.jld2")["times2"][1:T]
nve_p=load("PACF_JLD2/TIP3P_JLD2/vr_1_pacf.jld2")["pacf_mean"][1:T]
sd_001_p=load("PACF_JLD2/TIP3P_JLD2/sd_001_pacf.jld2")["pacf_mean"][1:T]
sd_01_p=load("PACF_JLD2/TIP3P_JLD2/sd_01_pacf.jld2")["pacf_mean"][1:T]
sd_1_p=load("PACF_JLD2/TIP3P_JLD2/sd_1_pacf.jld2")["pacf_mean"][1:T]
sd_10_p=load("PACF_JLD2/TIP3P_JLD2/sd_10_pacf.jld2")["pacf_mean"][1:T]
sd_100_p=load("PACF_JLD2/TIP3P_JLD2/sd_100_pacf.jld2")["pacf_mean"][1:T]
nve_p_int=integrate(time_p,abs.(nve_p))

p_001=integrate(time_p,abs.(sd_001_p.-nve_p))/nve_p_int
p_01=integrate(time_p,abs.(sd_01_p.-nve_p))/nve_p_int
p_1=integrate(time_p,abs.(sd_1_p.-nve_p))/nve_p_int
p_10=integrate(time_p,abs.(sd_10_p.-nve_p))/nve_p_int
p_100=integrate(time_p,abs.(sd_100_p.-nve_p))/nve_p_int
tip3p_sd_p=[p_001,p_01,p_1,p_10,p_100]


tau=[0.01,0.1,1,10,100]
tau2=[0.001,0.01, 0.1 , 1] 


T=200
T2=200

time_p=load("PACF_JLD2/TIP3P_JLD2/vr_1_pacf.jld2")["times2"][1:T]
nve_p=load("PACF_JLD2/TIP3P_JLD2/vr_1_pacf.jld2")["pacf_mean"][1:T]
sd_001_p=load("PACF_JLD2/TIP3P_JLD2/sd_001_pacf.jld2")["pacf_mean"][1:T]
sd_01_p=load("PACF_JLD2/TIP3P_JLD2/sd_01_pacf.jld2")["pacf_mean"][1:T]
sd_1_p=load("PACF_JLD2/TIP3P_JLD2/sd_1_pacf.jld2")["pacf_mean"][1:T]
sd_10_p=load("PACF_JLD2/TIP3P_JLD2/sd_10_pacf.jld2")["pacf_mean"][1:T]
sd_100_p=load("PACF_JLD2/TIP3P_JLD2/sd_100_pacf.jld2")["pacf_mean"][1:T]
nve_p_int=integrate(time_p,abs.(nve_p))

p_001=integrate(time_p,abs.(sd_001_p.-nve_p))/nve_p_int
p_01=integrate(time_p,abs.(sd_01_p.-nve_p))/nve_p_int
p_1=integrate(time_p,abs.(sd_1_p.-nve_p))/nve_p_int
p_10=integrate(time_p,abs.(sd_10_p.-nve_p))/nve_p_int
p_100=integrate(time_p,abs.(sd_100_p.-nve_p))/nve_p_int
sd_tip3p=[p_001,p_01,p_1,p_10,p_100]


time_p=load("PACF_JLD2/TIP4P_JLD2/vr_1_pacf.jld2")["times2"][1:T]
nve_p=load("PACF_JLD2/TIP4P_JLD2/vr_1_pacf.jld2")["pacf_mean"][1:T]
sd_001_p=load("PACF_JLD2/TIP4P_JLD2/sd_001_pacf.jld2")["pacf_mean"][1:T]
sd_01_p=load("PACF_JLD2/TIP4P_JLD2/sd_01_pacf.jld2")["pacf_mean"][1:T]
sd_1_p=load("PACF_JLD2/TIP4P_JLD2/sd_1_pacf.jld2")["pacf_mean"][1:T]
sd_10_p=load("PACF_JLD2/TIP4P_JLD2/sd_10_pacf.jld2")["pacf_mean"][1:T]
sd_100_p=load("PACF_JLD2/TIP4P_JLD2/sd_100_pacf.jld2")["pacf_mean"][1:T]
nve_p_int=integrate(time_p,abs.(nve_p))

p_001=integrate(time_p,abs.(sd_001_p.-nve_p))/nve_p_int
p_01=integrate(time_p,abs.(sd_01_p.-nve_p))/nve_p_int
p_1=integrate(time_p,abs.(sd_1_p.-nve_p))/nve_p_int
p_10=integrate(time_p,abs.(sd_10_p.-nve_p))/nve_p_int
p_100=integrate(time_p,abs.(sd_100_p.-nve_p))/nve_p_int
sd_tip4p=[p_001,p_01,p_1,p_10,p_100]


time_p=load("PACF_JLD2/Anilin_JLD2/vr_1_pacf.jld2")["times2"][1:T]
nve_p=load("PACF_JLD2/Anilin_JLD2/vr_1_pacf.jld2")["pacf_mean"][1:T]
sd_001_p=load("PACF_JLD2/Anilin_JLD2/sd_001_pacf.jld2")["pacf_mean"][1:T]
sd_01_p=load("PACF_JLD2/Anilin_JLD2/sd_01_pacf.jld2")["pacf_mean"][1:T]
sd_1_p=load("PACF_JLD2/Anilin_JLD2/sd_1_pacf.jld2")["pacf_mean"][1:T]
sd_10_p=load("PACF_JLD2/Anilin_JLD2/sd_10_pacf.jld2")["pacf_mean"][1:T]
sd_100_p=load("PACF_JLD2/Anilin_JLD2/sd_100_pacf.jld2")["pacf_mean"][1:T]
nve_p_int=integrate(time_p,abs.(nve_p))

p_001=integrate(time_p,abs.(sd_001_p.-nve_p))/nve_p_int
p_01=integrate(time_p,abs.(sd_01_p.-nve_p))/nve_p_int
p_1=integrate(time_p,abs.(sd_1_p.-nve_p))/nve_p_int
p_10=integrate(time_p,abs.(sd_10_p.-nve_p))/nve_p_int
p_100=integrate(time_p,abs.(sd_100_p.-nve_p))/nve_p_int
sd_anilin=[p_001,p_01,p_1,p_10,p_100]

time_p=load("PACF_JLD2/PENTANE_JLD2/vr_1_pacf.jld2")["times2"][1:T]
nve_p=load("PACF_JLD2/PENTANE_JLD2/vr_1_pacf.jld2")["pacf_mean"][1:T]
sd_001_p=load("PACF_JLD2/PENTANE_JLD2/sd_001_pacf.jld2")["pacf_mean"][1:T]
sd_01_p=load("PACF_JLD2/PENTANE_JLD2/sd_01_pacf.jld2")["pacf_mean"][1:T]
sd_1_p=load("PACF_JLD2/PENTANE_JLD2/sd_1_pacf.jld2")["pacf_mean"][1:T]
sd_10_p=load("PACF_JLD2/PENTANE_JLD2/sd_10_pacf.jld2")["pacf_mean"][1:T]
sd_100_p=load("PACF_JLD2/PENTANE_JLD2/sd_100_pacf.jld2")["pacf_mean"][1:T]
nve_p_int=integrate(time_p,abs.(nve_p))

p_001=integrate(time_p,abs.(sd_001_p.-nve_p))/nve_p_int
p_01=integrate(time_p,abs.(sd_01_p.-nve_p))/nve_p_int
p_1=integrate(time_p,abs.(sd_1_p.-nve_p))/nve_p_int
p_10=integrate(time_p,abs.(sd_10_p.-nve_p))/nve_p_int
p_100=integrate(time_p,abs.(sd_100_p.-nve_p))/nve_p_int
sd_pentane=[p_001,p_01,p_1,p_10,p_100]

time_p=load("PACF_JLD2/GLYCERINE_JLD2/vr_1_pacf.jld2")["times2"][1:T2]
nve_p=load("PACF_JLD2/GLYCERINE_JLD2/vr_1_pacf.jld2")["pacf_mean"][1:T2]
sd_001_p=load("PACF_JLD2/GLYCERINE_JLD2/sd_001_pacf.jld2")["pacf_mean"][1:T2]
sd_01_p=load("PACF_JLD2/GLYCERINE_JLD2/sd_01_pacf.jld2")["pacf_mean"][1:T2]
sd_1_p=load("PACF_JLD2/GLYCERINE_JLD2/sd_1_pacf.jld2")["pacf_mean"][1:T2]
sd_10_p=load("PACF_JLD2/GLYCERINE_JLD2/sd_10_pacf.jld2")["pacf_mean"][1:T2]
sd_100_p=load("PACF_JLD2/GLYCERINE_JLD2/sd_100_pacf.jld2")["pacf_mean"][1:T2]
nve_p_int=integrate(time_p,abs.(nve_p))

p_001=integrate(time_p,abs.(sd_001_p.-nve_p))/nve_p_int
p_01=integrate(time_p,abs.(sd_01_p.-nve_p))/nve_p_int
p_1=integrate(time_p,abs.(sd_1_p.-nve_p))/nve_p_int
p_10=integrate(time_p,abs.(sd_10_p.-nve_p))/nve_p_int
p_100=integrate(time_p,abs.(sd_100_p.-nve_p))/nve_p_int
sd_glycerine=[p_001,p_01,p_1,p_10,p_100]


vr_1=[0.95332916, 0.93254085, 0.91324553, 0.89473093, 0.876904  ,
0.85975644, 0.84321841, 0.82727585, 0.81189535, 0.79709845,
0.78287278, 0.7692063 , 0.75605103, 0.74339258, 0.73119415,
0.71940745, 0.70804298, 0.69714244, 0.68667616, 0.67657991,
0.6668347 , 0.65740205, 0.64833892, 0.63954919, 0.63104084,
0.62281284, 0.61482131, 0.60711419, 0.59966133, 0.59247442,
0.58554974, 0.57889506, 0.57245617, 0.56629051, 0.56025029,
0.55436685, 0.5486614 , 0.54320375, 0.53795121, 0.53285647,
0.5279036 , 0.52311854, 0.51851638, 0.51409626, 0.50981479,
0.50569232, 0.50168545, 0.49782825, 0.49414394, 0.49057388,
0.48711518, 0.48376987, 0.48054576, 0.47743454, 0.4744013 ,
0.4714389 , 0.46857041, 0.46580751, 0.46309615, 0.46051463]

sd_001=[0.94410675, 0.93885412, 0.93393054, 0.92905727, 0.92424495,
0.91950067, 0.91482139, 0.91020136, 0.90561662, 0.90107477,
0.89656685, 0.89209759, 0.8876702 , 0.88328536, 0.8789385 ,
0.8746271 , 0.87035925, 0.86613513, 0.86195795, 0.85781572,
0.85371636, 0.84966513, 0.84564901, 0.84167847, 0.83773562,
0.83382063, 0.82994413, 0.82610444, 0.82229868, 0.81853344,
0.81479488, 0.8110808 , 0.80740757, 0.8037752 , 0.80018048,
0.79660636, 0.7930665 , 0.78955096, 0.78608253, 0.7826438 ,
0.77922889, 0.7758452 , 0.7724963 , 0.76917964, 0.76589185,
0.76264779, 0.75943682, 0.75625809, 0.75310874, 0.74997864,
0.74688113, 0.74379711, 0.74073302, 0.73769358, 0.73468943,
0.73171448, 0.72877127, 0.72585964, 0.72296927, 0.72010508]

sd_01=[0.98314223, 0.96283096, 0.94337097, 0.92460243, 0.90649894,
0.88905507, 0.87225128, 0.85609603, 0.84053602, 0.82549178,
0.81095716, 0.79697434, 0.78352066, 0.77052119, 0.75803432,
0.74603122, 0.73446492, 0.72328567, 0.71242382, 0.70194058,
0.69178565, 0.68194602, 0.67244894, 0.66331534, 0.65454047,
0.64606771, 0.63783727, 0.62992828, 0.62230012, 0.61492273,
0.60781508, 0.60092585, 0.59430319, 0.58792075, 0.58176006,
0.57579455, 0.57003092, 0.56444508, 0.55897406, 0.5536912 ,
0.54860301, 0.54372214, 0.53893518, 0.53428591, 0.52982989,
0.52553073, 0.52136206, 0.51730522, 0.51340119, 0.50964295,
0.50603593, 0.50256063, 0.49926662, 0.49609307, 0.49303207,
0.4900764 , 0.48723169, 0.48450304, 0.48181782, 0.47921455]

sd_1=[0.97816079, 0.95492437, 0.93318716, 0.91239249, 0.89243146,
0.8733057 , 0.8550016 , 0.83743621, 0.82056852, 0.80438037,
0.78880751, 0.77383197, 0.75939676, 0.74551985, 0.73216603,
0.71927287, 0.70687466, 0.69492566, 0.68339176, 0.67230871,
0.66164403, 0.65139734, 0.64152926, 0.63200658, 0.62282713,
0.61394497, 0.60531964, 0.59698562, 0.58893928, 0.58122581,
0.57378968, 0.56663596, 0.55975559, 0.55313277, 0.54681287,
0.54071278, 0.53477658, 0.5290727 , 0.52356303, 0.51827406,
0.51320689, 0.5083372 , 0.50365245, 0.49919495, 0.49490025,
0.49077217, 0.48674537, 0.48287085, 0.47914571, 0.47554471,
0.4720586 , 0.46865072, 0.46538841, 0.46226676, 0.45927452,
0.4564008 , 0.45362291, 0.45096881, 0.44841   , 0.44593522]

sd_10=[0.97057518, 0.94842735, 0.92790569, 0.90825234, 0.88935221,
0.87117783, 0.85370674, 0.83691126, 0.82075626, 0.80515985,
0.79012488, 0.77565885, 0.76171609, 0.74828088, 0.73533772,
0.72289267, 0.71095215, 0.6994812 , 0.68849265, 0.67791033,
0.66770445, 0.65785342, 0.64831942, 0.63911921, 0.63028331,
0.62178852, 0.61356615, 0.60558892, 0.59787269, 0.59038054,
0.58314476, 0.57614287, 0.56932438, 0.56277564, 0.55647256,
0.55042782, 0.54462463, 0.53902858, 0.53367801, 0.52850906,
0.52346803, 0.51860201, 0.51390156, 0.50935007, 0.50493096,
0.50063084, 0.49648451, 0.49246396, 0.48857739, 0.48481945,
0.48114911, 0.47760187, 0.47415418, 0.47085189, 0.46763987,
0.46453864, 0.46152466, 0.45861967, 0.45581746, 0.45310301]

sd_100=[0.96544269, 0.94466915, 0.92544712, 0.90705161, 0.88937278,
0.87235606, 0.85594948, 0.84019357, 0.82503482, 0.81043753,
0.79632963, 0.78272793, 0.76967054, 0.75706599, 0.74496158,
0.73325388, 0.72191257, 0.71099674, 0.70049425, 0.69031605,
0.68047879, 0.67099441, 0.6618158 , 0.65292874, 0.64428405,
0.63594442, 0.62790103, 0.62008353, 0.61254998, 0.60530469,
0.5983293 , 0.59159004, 0.58504618, 0.57873326, 0.57266182,
0.56680433, 0.56116824, 0.55576619, 0.55061568, 0.54560134,
0.54074066, 0.53606378, 0.53153242, 0.52714121, 0.52292393,
0.51880349, 0.51479686, 0.51089106, 0.50712489, 0.50346318,
0.4998744 , 0.496401  , 0.49300712, 0.48970747, 0.48649791,
0.48336371, 0.4803745 , 0.47751956, 0.47472109, 0.47198828]

lag=collect(1:60)
m_vr_1_int=integrate(lag,abs.(vr_1))
m_001=integrate(lag,abs.(sd_001.-vr_1))/m_vr_1_int
m_01=integrate(lag,abs.(sd_01.-vr_1))/m_vr_1_int
m_1=integrate(lag,abs.(sd_1.-vr_1))/m_vr_1_int
m_10=integrate(lag,abs.(sd_10.-vr_1))/m_vr_1_int
m_100=integrate(lag,abs.(sd_100.-vr_1))/m_vr_1_int

msn=[m_001,m_01,m_1,m_10,m_100]
sd_c= [NaN,1.4805665957762753, 0.679143411287748, 0.1189542494877709, 0.027543627563339227]

color1=["#D32F2F","#FF8C00","#D32F2F","#800080","#00509E"]
color2=["#D32F2F","#800080","#00509E","#006400","#FFCC00"]
TT=[9,18,27,35,43]
TT2=[0,4,9,18,27]

fig, ax = plt.subplots(figsize=(6, 3))
plt.bar(TT2.-3,sd_c,width=1,color=color1,alpha=0.95,edgecolor="black",linewidth=0.2, zorder=3)
plt.bar(TT.-2,tip3p_sd_v,width=1,color=color2,alpha=0.1,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar(TT.-1,sd_pentane,width=1,color=color2,alpha=0.25,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar(TT,sd_tip3p,width=1,color=color2,alpha=0.4,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar(TT.+1,sd_anilin,width=1,color=color2,alpha=0.55,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar(TT.+2,sd_glycerine,width=1,color=color2,alpha=0.7,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar(TT.+3,msn,width=1,color=color2,alpha=0.85,linewidth=0.2,edgecolor="black", zorder=3)

plt.bar(TT2.-3,sd_c,width=1,color="None",edgecolor="black",linewidth=0.2, zorder=4)
plt.bar(TT.-2,tip3p_sd_v,width=1,color="None",linewidth=0.2,edgecolor="black", zorder=4)
plt.bar(TT.-1,sd_pentane,width=1,color="None",linewidth=0.2,edgecolor="black", zorder=4)
plt.bar(TT,sd_tip3p,width=1,color="None",linewidth=0.2,edgecolor="black", zorder=4)
plt.bar(TT.+1,sd_anilin,width=1,color="None",linewidth=0.2,edgecolor="black", zorder=4)
plt.bar(TT.+2,sd_glycerine,width=1,color="None",linewidth=0.2,edgecolor="black", zorder=4)
plt.bar(TT.+3,msn,width=1,color="None",linewidth=0.2,edgecolor="black", zorder=4)

plt.bar([NaN,NaN],[NaN,NaN],label="vDOS",width=1,color="grey",alpha=0.95,edgecolor="black",linewidth=0.2, zorder=3)
plt.bar([NaN,NaN],[NaN,NaN],label="D(tip3p)",width=1,color="grey",alpha=0.1,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar([NaN,NaN],[NaN,NaN],label="η(pentane)",width=1,color="grey",alpha=0.25,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar([NaN,NaN],[NaN,NaN],label="η(tip3p)",width=1,color="grey",alpha=0.4,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar([NaN,NaN],[NaN,NaN],label="η(aninline)",width=1,color="grey",alpha=0.55,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar([NaN,NaN],[NaN,NaN],label="η(glycerol)",width=1,color="grey",alpha=0.7,linewidth=0.2,edgecolor="black", zorder=3)
plt.bar([NaN,NaN],[NaN,NaN],label="MSM",width=1,color="grey",alpha=0.85,linewidth=0.2,edgecolor="black", zorder=3)

ax.grid(axis="y",zorder=0,linestyle="dashed",linewidth=0.5)
ax.set_yscale("log")
ax.set_xlabel(L"$\tau_T$/lag time in ps", fontsize=14)
ax.set_ylabel("Relative Error",fontsize=14)
ax.set_xticks([1,9,18,27,35,43])
ax.set_yticks([0.01,0.1,1])
ax.set_ylim([0.005,10])
ax.set_xlim([0,47])
ax.set_yticklabels(["1%","10%","100%"])
ax.set_xticklabels(["0.001 ps","0.01 ps","0.1 ps","1 ps","10 ps","100 ps"])
ax.tick_params(axis="both",which="major",labelsize=12)
plt.tight_layout()
plt.legend(fontsize=12,ncol=2,borderpad=0.3,handletextpad=0.3,handlelength=1,columnspacing=0.7,framealpha=1)
plt.savefig("error_plot.pdf")
#plt.show()
plt.close()

println("End")