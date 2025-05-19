using Chemfiles
using FFTW
using StatsBase
using NumericalIntegration
using JLD2


println("Start")

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
function VAC(Values)
    N=length(Values)
    X=zeros(2*N)
    X[1:N]=Values
    X=fft(X)
    X=X.*conj.(X)
    X=ifft(X)
    X=real.(X[1:N])./N
    return X
end

function create_vacf(name,outname)
    traj=Trajectory(name)
    frames=Int64(length(traj))
    times=collect(1:10:10*frames)*1E-12
    atoms=floor(Int64,length(read_step(traj,0))/3)
    vel=zeros(frames,3,atoms)*1000
    i=0
    while i < frames
        vel_atom=velocities(read_step(traj,i))
        vel[i+1,:,:]=center_of_mass(vel_atom)
        i+=1
    end

    VACF=zeros(frames)
    for i in collect(1:atoms)
        try
            VACF+=VAC(vel[:,1,i])/3+VAC(vel[:,2,i])/3+VAC(vel[:,3,i])/3
        catch
            println(i)
        end
    end
    VACF=VACF/atoms
    jldopen("$outname.jld2", "w+") do file
        file["VACF"]=VACF
        file["VACF_N"]=VACF./VACF[1]
        file["times"]=times
        file["VACF_C"]=cumul_integrate(times,VACF)
        file["VACF_CN"]=cumul_integrate(times,VACF./VACF[1])
    end
    println("Finished with: $outname")
    GC.gc()
end

#create_vacf("VR_100/vr_100_vacf.trr","vr_100")
#create_vacf("VR_10/vr_10_vacf.trr","vr_10")
#create_vacf("VR_1/vr_1_vacf.trr","vr_1")
#create_vacf("Velocity_Rescale/velocityrescale_vacf.trr","vr_01")
#create_vacf("VR_01/vr_01_vacf.trr","vr_001")

#create_vacf("SD_2/sd_2_vacf.trr","sd_2")
#create_vacf("Nose_Hover/nose_hover_vacf.trr","nh")
#create_vacf("NVE/nve_vacf2.trr","nve")
create_vacf("Berendsen/berendsen_vacf.trr","be")
