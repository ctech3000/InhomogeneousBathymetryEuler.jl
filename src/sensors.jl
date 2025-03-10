mutable struct Sensors
    sensors_pos_x::Vector{Float64}
    sensors_pos_χ::Vector{Float64}
    data::Vector{Vector{Float64}}
end

function Sensors(pos::Vector{Float64},domain::AbstractDomain,trans::σTransform,time_vec::Vector{Float64},var::String)
    nSensors = length(pos)
    nTime = length(time_vec)
    if var == "x"
        for x in pos
            if x < domain.x_L || x > domain.x_R
                print("Error in Sensors: Position out of bounds!\n")
                break
            end
        end
        sensors_pos_x = pos
        sensors_pos_χ = trans.χ(pos)
    elseif var == "χ"
        χ_L = trans.χ(domain.x_L)
        χ_R = trans.χ(domain.x_R)
        for χ in pos
            if χ < χ_R || χ > χ_L
                print("Error in Sensors: Position out of bounds!\n")
                break
            end
        end
        sensors_pos_x = trans.x(pos)
        sensors_pos_χ = pos
    else
        print("In Sensors: Invalid var!\n")
    end
    data = [zeros(Float64,nTime) for s = 1:nSensors]
    
    return Sensors(sensors_pos_x,sensors_pos_χ,data)
end

function Sensors(domain::AbstractDomain,trans::σTransform,time_vec::Vector{Float64})
    pos = [0.0,2.0,4.0,6.0]
    return Sensors(pos,domain,trans,time_vec,"x")
end

function extractSensorData!(sensors::Sensors,etas::Vector{Vector{Float64}},χs::Vector{Float64}; save_phi::Union{Bool,Tuple{Int64, Int64, Int64}}=false)
    nSensors = length(sensors.sensors_pos_x)
    interp_points = Vector{Any}(undef,nSensors)
    if isa(save_phi,Tuple{Int64, Int64, Int64})
        loc_χs = χs[1:save_phi[1]:end]
        oldNTime = length(sensors.data[1])
        newNTime = round(Integer,(oldNTime-1)/save_phi[3]+1)
        sensors.data = [zeros(Float64,newNTime) for s = 1:nSensors]
    else
        loc_χs = χs
    end
    for (sensor_idx,sensor_pos) in enumerate(sensors.sensors_pos_χ)
        i = 1
        while sensor_pos > loc_χs[i]
            i += 1
        end
        if sensor_pos == loc_χs[i]
            interp_points[sensor_idx] = [(i,1.0),(i,0.0)]
        else
            d1 = sensor_pos - loc_χs[i-1]
            d2 = loc_χs[i] - sensor_pos
            interp_points[sensor_idx] = [(i-1,d1),(i,d2)]
        end
    end
    for (t_idx,eta) in enumerate(etas)
        for sensor_idx = 1:nSensors
            sensors.data[sensor_idx][t_idx] = (interp_points[sensor_idx][2][2]*eta[interp_points[sensor_idx][1][1]] + interp_points[sensor_idx][1][2]*eta[interp_points[sensor_idx][2][1]])/(interp_points[sensor_idx][1][2]+interp_points[sensor_idx][2][2])
        end
    end
end