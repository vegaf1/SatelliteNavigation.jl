module SatelliteNavigation

#using GLMakie, FileIO, Colors
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using SatelliteDynamics
using DelimitedFiles
using ForwardDiff
using PlotlyJS
using Distributions

include("Doppler4.jl")

# function test()
#     println("test function works")
# end

function create_sat_orbit(orbit_params)

    #orbit_params[1] -> true anomaly seperation
    #orbit_params[2] -> RAAN seperation
    #orbit_params[3] -> Delta True Anomaly seperation
    #orbit_params[4] -> altitude
    #orbit_params[5] -> Number of Satellites

    iss1 = [6371e3 + (orbit_params[4]*1e3), 0.0004879, 90.6391, 194.5859- (orbit_params[2]/2), 151.2014, 190];
    iss2 = [6371e3 + (orbit_params[4]*1e3), 0.0004879, 90.6391, 194.5859 - (orbit_params[2]/2), 151.2014, 190+orbit_params[1]];
    iss3 = [6371e3 + (orbit_params[4]*1e3), 0.0004879, 90.6391, 194.5859 + (orbit_params[2]/2), 151.2014, 190+orbit_params[3]]; 
    iss4 = [6371e3 + (orbit_params[4]*1e3), 0.0004879, 90.6391, 194.5859 + (orbit_params[2]/2), 151.2014, 190+orbit_params[3]+orbit_params[1]]; 
    

    eci0_1 = sOSCtoCART(iss1, use_degrees=true)
    eci0_2 = sOSCtoCART(iss2, use_degrees=true)
    eci0_3 = sOSCtoCART(iss3, use_degrees=true)
    eci0_4 = sOSCtoCART(iss4, use_degrees=true)

    #orbital period is the same since all have the same altitude
    T = orbit_period(iss1[1])

    # Declare simulation initial Epoch. Start time
    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)  #year, month, day, hour, minute, seconds, nanoseconds 


    #final time for one orbit. Adds initial time to the orbital period
    epcf = epc0 + T

    # Create an EarthInertialState orbit propagagator
    #needs initial epoch of state and the state vector
    orb1  = EarthInertialState(epc0, eci0_1, dt=1.0,
                               mass=1.0, n_grav=0, m_grav=0,
                               drag=false, srp=false,
                               moon=false, sun=false,
                               relativity=false)

    orb2  = EarthInertialState(epc0, eci0_2, dt=1.0,
                               mass=1.0, n_grav=0, m_grav=0,
                               drag=false, srp=false,
                               moon=false, sun=false,
                               relativity=false)

    orb3  = EarthInertialState(epc0, eci0_3, dt=1.0,
                               mass=1.0, n_grav=0, m_grav=0,
                               drag=false, srp=false,
                               moon=false, sun=false,
                               relativity=false)

    orb4  = EarthInertialState(epc0, eci0_4, dt=1.0,
                               mass=1.0, n_grav=0, m_grav=0,
                               drag=false, srp=false,
                               moon=false, sun=false,
                               relativity=false)

    # Propagate the orbit
    # orbit until the final time

    t_1, epc_1, eci_1 = sim!(orb1, epcf);
    t_2, epc_2, eci_2 = sim!(orb2, epcf);
    t_3, epc_3, eci_3 = sim!(orb3, epcf);
    t_4, epc_4, eci_4 = sim!(orb4, epcf);

    #add 3 additional satellite trajectories if the number of satellites is 7
    if orbit_params[5] ==7

        iss5 = [6371e3 + (orbit_params[4]*1e3), 0.0004879, 90.6391, 194.5859+(1.5*orbit_params[2]), 151.2014, 190+(2*orbit_params[3])]; 
        iss6 = [6371e3 + (orbit_params[4]*1e3), 0.0004879, 90.6391, 194.5859+(1.5*orbit_params[2]), 151.2014, 190+(2*orbit_params[3])+orbit_params[1]]; 
        iss7 = [6371e3 + (orbit_params[4]*1e3), 0.0004879, 90.6391, 194.5859+(1.5*orbit_params[2]), 151.2014, 190+(2*orbit_params[3])+(2*orbit_params[1])];

        eci0_5 = sOSCtoCART(iss5, use_degrees=true)
        eci0_6 = sOSCtoCART(iss6, use_degrees=true)
        eci0_7 = sOSCtoCART(iss7, use_degrees=true)

        orb5  = EarthInertialState(epc0, eci0_5, dt=1.0,
                                   mass=1.0, n_grav=0, m_grav=0,
                                   drag=false, srp=false,
                                   moon=false, sun=false,
                                   relativity=false)

        orb6  = EarthInertialState(epc0, eci0_6, dt=1.0,
                                   mass=1.0, n_grav=0, m_grav=0,
                                   drag=false, srp=false,
                                   moon=false, sun=false,
                                   relativity=false)

        orb7  = EarthInertialState(epc0, eci0_7, dt=1.0,
                                   mass=1.0, n_grav=0, m_grav=0,
                                   drag=false, srp=false,
                                   moon=false, sun=false,
                                   relativity=false)

        #simulate additional 3 satellites
        t_5, epc_5, eci_5 = sim!(orb5, epcf);
        t_6, epc_6, eci_6 = sim!(orb6, epcf);
        t_7, epc_7, eci_7 = sim!(orb7, epcf);

    end

    #return the poses of each satellite depending on the user selecting 4 or 7 satellites
    if orbit_params[5] == 7

        combined_eci = [eci_1; eci_2; eci_3; eci_4; eci_5; eci_6; eci_7]
    else
        combined_eci = [eci_1; eci_2; eci_3; eci_4]
    end

    return combined_eci

end

function plot_sat_constellation(combined_eci, tag, satnum)

    #Plot Satellite 1 and 2 Orbit
    sat1 = scatter(x=combined_eci[1,:], y=combined_eci[2,:], z=combined_eci[3,:], type="scatter3d", mode="lines", name="orbit 1&2")

    #Plot Satellite 3 and 4 Orbit
    sat3 = scatter(x=combined_eci[13,:], y=combined_eci[14,:], z=combined_eci[15,:], type="scatter3d", mode="lines", name="orbit 3&4")

    #Plot all satellite initial positions
    satellites_4 = scatter(x=[combined_eci[1,1],combined_eci[7,1], combined_eci[13,1], combined_eci[19,1]], y=[combined_eci[2,1],combined_eci[8,1], combined_eci[14,1], combined_eci[20,1]],z=[combined_eci[3,1],combined_eci[9,1], combined_eci[15,1], combined_eci[21,1]], mode="markers", marker_size = 4, type="scatter3d", name="satellites")

    #Plot Tag position
    tag = scatter(x = [tag[1]], y = [tag[2]], z = [tag[3]], type="scatter3d", name="tag", mode="markers", marker_size=5)

    #Plot the Guess
    #guess = scatter(x = [xyz[1], xyz[1]], y = [xyz[2], xyz[2]], z = [xyz[3], xyz[3]], type="scatter3d", name="guess", mode="markers", marker_size=5)

    #Plotting the sphere (Earth)
    n = 100
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    x = cos.(u) * sin.(v)'
    y = sin.(u) * sin.(v)' 
    z = ones(n) * cos.(v)' 

    earth = surface(z=z*6371000, x=x*6371000, y=y*6371000, showscale = false)


    #add extra orbit and satellites to the plot if there is 7
    if satnum == 7

        sat5 = scatter(x=combined_eci[25,:], y=combined_eci[26,:], z=combined_eci[27,:], type="scatter3d", mode="lines", name="orbit 5&6&7")
        satellites_7 = scatter(x=[combined_eci[1,1],combined_eci[7,1], combined_eci[13,1], combined_eci[19,1], combined_eci[25,1], combined_eci[31,1], combined_eci[37,1]], y=[combined_eci[2,1],combined_eci[8,1], combined_eci[14,1], combined_eci[20,1], combined_eci[26,1], combined_eci[32,1], combined_eci[38,1]],z=[combined_eci[3,1],combined_eci[9,1], combined_eci[15,1], combined_eci[21,1], combined_eci[27,1], combined_eci[33,1], combined_eci[39,1]], mode="markers", marker_size = 4, type="scatter3d", name="satellites")
        plot([sat1, sat3, sat5, satellites_7, tag, earth])

    else

        plot([sat1, sat3, satellites_4, tag, earth])

    end

end

#find the zenith angle for each satellite
function zenith(satposes, r0, satnum)

    zenithangle = zeros(satnum)


    #vector between tag and the satellite
    normalvec = r0
  
    for i in 1:satnum

        vector = satposes[i,:] - r0[1:3]
        
        #Find the angle between the normal and the vector going from the tag to the satellte
        theta = acos(dot(normalvec,vector)/(norm(normalvec)*norm(vector))) * (180/pi)
        
        #save elevation angle
        zenithangle[i] = theta
    end

    return zenithangle

end

#FUNCTIONS FOR TOA ANALYSIS

function measurment_residual(r0, satposes, zenith, t, frequency)

    hi = 350e3 #in meters
    distance_scale = R_EARTH*(1e-3) # scales km to custom scale
    time_scale = 1/(C_LIGHT/R_EARTH) #scales seconds to custom scale
    c = 1 # speed of light

    omega = OMEGA_EARTH*time_scale
    
    scale_f1 = (40.3/(frequency[1])^2)*1e22
    scale_f2 = (40.3/(frequency[2])^2)*1e22
    
    angles = zeros(4)
    
    for i=1:4
        angles[i] = asind((R_EARTH*sind(zenith[i]))/(R_EARTH + hi))
    end
    
    A1 = [cos(omega*t[1]) sin(omega*t[1]) 0;
          -sin(omega*t[1]) cos(omega*t[1]) 0;
          0 0 1];
    A2 = [cos(omega*t[2]) sin(omega*t[2]) 0;
          -sin(omega*t[2]) cos(omega*t[2]) 0;
          0 0 1];
    A3 = [cos(omega*t[3]) sin(omega*t[3]) 0;
          -sin(omega*t[3]) cos(omega*t[3]) 0;
          0 0 1];
    A4 = [cos(omega*t[4]) sin(omega*t[4]) 0;
          -sin(omega*t[4]) cos(omega*t[4]) 0;
          0 0 1];

    #if there are two frequency, calculate additional matrices
    if frequency[2] != 0

        A5 = [cos(omega*t[5]) sin(omega*t[5]) 0;
            -sin(omega*t[5]) cos(omega*t[5]) 0;
            0 0 1];
        A6 = [cos(omega*t[6]) sin(omega*t[6]) 0;
            -sin(omega*t[6]) cos(omega*t[6]) 0;
            0 0 1];
        A7 = [cos(omega*t[7]) sin(omega*t[7]) 0;
            -sin(omega*t[7]) cos(omega*t[7]) 0;
            0 0 1];
        A8 = [cos(omega*t[8]) sin(omega*t[8]) 0;
            -sin(omega*t[8]) cos(omega*t[8]) 0;
            0 0 1];
    end
    
    #for 1 frequency case, we insert the noise into the time measurment. Direct effect on pseudo range
    res1 = norm(r0[1:3] - A1*satposes[1,1:3]) - c*(t[1] - r0[4]);
    res2 = norm(r0[1:3] - A2*satposes[2,1:3]) - c*(t[2] - r0[4]);
    res3 = norm(r0[1:3] - A3*satposes[3,1:3]) - c*(t[3] - r0[4]);
    res4 = norm(r0[1:3] - A4*satposes[4,1:3]) - c*(t[4] - r0[4]);
    
    residuals = [res1, res2, res3, res4]

    if frequency[2] != 0

        res1 = norm(r0[1:3] - A1*satposes[1,1:3]) - c*(t[1] - r0[4]) + (scale_f1*(r0[5]/cosd(angles[1])))*(1e-3)/distance_scale;
        res2 = norm(r0[1:3] - A2*satposes[2,1:3]) - c*(t[2] - r0[4]) + (scale_f1*(r0[6]/cosd(angles[2])))*(1e-3)/distance_scale;
        res3 = norm(r0[1:3] - A3*satposes[3,1:3]) - c*(t[3] - r0[4]) + (scale_f1*(r0[7]/cosd(angles[3])))*(1e-3)/distance_scale;
        res4 = norm(r0[1:3] - A4*satposes[4,1:3]) - c*(t[4] - r0[4]) + (scale_f1*(r0[8]/cosd(angles[4])))*(1e-3)/distance_scale;
        res5 = norm(r0[1:3] - A5*satposes[1,1:3]) - c*(t[5] - r0[4]) + (scale_f2*(r0[5]/cosd(angles[1])))*(1e-3)/distance_scale;
        res6 = norm(r0[1:3] - A6*satposes[2,1:3]) - c*(t[6] - r0[4]) + (scale_f2*(r0[6]/cosd(angles[2])))*(1e-3)/distance_scale;
        res7 = norm(r0[1:3] - A7*satposes[3,1:3]) - c*(t[7] - r0[4]) + (scale_f2*(r0[7]/cosd(angles[3])))*(1e-3)/distance_scale;
        res8 = norm(r0[1:3] - A8*satposes[4,1:3]) - c*(t[8] - r0[4]) + (scale_f2*(r0[8]/cosd(angles[4])))*(1e-3)/distance_scale;

        residuals = [res1, res2, res3, res4, res5, res6, res7, res8]

    end

    return residuals
    
end

function TOA_residual(x, sat_pose_t, zenith, frequency)

    hi = 350e3 #in meters
    distance_scale = R_EARTH*(1e-3) # scales km to custom scale
    time_scale = 1/(C_LIGHT/R_EARTH) #scales seconds to custom scale
    c = 1 # speed of light

    omega = OMEGA_EARTH*time_scale
    
    angles = zeros(4)
    
    for i=1:4
        angles[i] = asind((R_EARTH*sind(zenith[i]))/(R_EARTH + hi))
    end

    scale_f1 = (40.3/(frequency[1])^2)*1e22
    scale_f2 = (40.3/(frequency[2])^2)*1e22
    
    A1 = [cos(omega*sat_pose_t[1,4]) sin(omega*sat_pose_t[1,4]) 0;
          -sin(omega*sat_pose_t[1,4]) cos(omega*sat_pose_t[1,4]) 0;
          0 0 1];
    A2 = [cos(omega*sat_pose_t[2,4]) sin(omega*sat_pose_t[2,4]) 0;
          -sin(omega*sat_pose_t[2,4]) cos(omega*sat_pose_t[2,4]) 0;
          0 0 1];
    A3 = [cos(omega*sat_pose_t[3,4]) sin(omega*sat_pose_t[3,4]) 0;
          -sin(omega*sat_pose_t[3,4]) cos(omega*sat_pose_t[3,4]) 0;
          0 0 1];
    A4 = [cos(omega*sat_pose_t[4,4]) sin(omega*sat_pose_t[4,4]) 0;
          -sin(omega*sat_pose_t[4,4]) cos(omega*sat_pose_t[4,4]) 0;
          0 0 1];

    if frequency[2] != 0

        A5 = [cos(omega*sat_pose_t[5,4]) sin(omega*sat_pose_t[5,4]) 0;
            -sin(omega*sat_pose_t[5,4]) cos(omega*sat_pose_t[5,4]) 0;
            0 0 1];
        A6 = [cos(omega*sat_pose_t[6,4]) sin(omega*sat_pose_t[6,4]) 0;
            -sin(omega*sat_pose_t[6,4]) cos(omega*sat_pose_t[6,4]) 0;
            0 0 1];
        A7 = [cos(omega*sat_pose_t[7,4]) sin(omega*sat_pose_t[7,4]) 0;
            -sin(omega*sat_pose_t[7,4]) cos(omega*sat_pose_t[7,4]) 0;
            0 0 1];
        A8 = [cos(omega*sat_pose_t[8,4]) sin(omega*sat_pose_t[8,4]) 0;
            -sin(omega*sat_pose_t[8,4]) cos(omega*sat_pose_t[8,4]) 0;
            0 0 1];
    end
    
    res1 = norm(x[1:3] - A1*sat_pose_t[1,1:3]) - c*(sat_pose_t[1,4] - x[4]);
    res2 = norm(x[1:3] - A2*sat_pose_t[2,1:3]) - c*(sat_pose_t[2,4] - x[4]);
    res3 = norm(x[1:3] - A3*sat_pose_t[3,1:3]) - c*(sat_pose_t[3,4] - x[4]);
    res4 = norm(x[1:3] - A4*sat_pose_t[4,1:3]) - c*(sat_pose_t[4,4] - x[4]);
    
    TOA_residual = [res1, res2, res3, res4]

    if frequency[2] != 0

        res1 = norm(x[1:3] - A1*sat_pose_t[1,1:3]) - c*(sat_pose_t[1,4] - x[4]) + (scale_f1*(x[5]/cosd(angles[1])))*(1e-3)/distance_scale;
        res2 = norm(x[1:3] - A2*sat_pose_t[2,1:3]) - c*(sat_pose_t[2,4] - x[4]) + (scale_f1*(x[6]/cosd(angles[2])))*(1e-3)/distance_scale;
        res3 = norm(x[1:3] - A3*sat_pose_t[3,1:3]) - c*(sat_pose_t[3,4] - x[4]) + (scale_f1*(x[7]/cosd(angles[3])))*(1e-3)/distance_scale;
        res4 = norm(x[1:3] - A4*sat_pose_t[4,1:3]) - c*(sat_pose_t[4,4] - x[4]) + (scale_f1*(x[8]/cosd(angles[4])))*(1e-3)/distance_scale;

        res5 = norm(x[1:3] - A5*sat_pose_t[1,1:3]) - c*(sat_pose_t[5,4] - x[4]) + (scale_f2*(x[5]/cosd(angles[1])))*(1e-3)/distance_scale;
        res6 = norm(x[1:3] - A6*sat_pose_t[2,1:3]) - c*(sat_pose_t[6,4] - x[4]) + (scale_f2*(x[6]/cosd(angles[2])))*(1e-3)/distance_scale;
        res7 = norm(x[1:3] - A7*sat_pose_t[3,1:3]) - c*(sat_pose_t[7,4] - x[4]) + (scale_f2*(x[7]/cosd(angles[3])))*(1e-3)/distance_scale;
        res8 = norm(x[1:3] - A8*sat_pose_t[4,1:3]) - c*(sat_pose_t[8,4] - x[4]) + (scale_f2*(x[8]/cosd(angles[4])))*(1e-3)/distance_scale;
        
        TOA_residual = [res1, res2, res3, res4, res5, res6, res7, res8]
    end

    return TOA_residual

end

function get_time(r0_TEC, satposes, zenith, frequency)

    distance_scale = R_EARTH*(1e-3) # scales km to custom scale
    time_scale = 1/(C_LIGHT/R_EARTH) #scales seconds to custom scale
    c = 1 # speed of light


    n = 1000 #max num of iterations

    time = [zeros(8) for i = 1:1000]
    Rt = [zeros(8) for i = 1:1000]

    #initial guess for time (backslash?)
    time[1] = [0.1\time_scale, 0.1\time_scale, 0.1\time_scale, 0.1\time_scale, 0.1\time_scale, 0.1\time_scale, 0.1\time_scale, 0.1\time_scale] #initial guess

    #get the measurment residual
    Rt[1] = measurment_residual(r0_TEC, satposes, zenith, time[1], frequency)

    iters = 0

    for i=1:n

        Rt[i] = measurment_residual(r0_TEC, satposes, zenith, time[i], frequency)

        #println("this is norm of residual: ", norm(Rt[i]))
        
        iters += 1

        if(norm(Rt[i]) < 1e-6)
            #println("time converged")
            break
        end

        jacobian = ForwardDiff.jacobian(dt -> measurment_residual(r0_TEC, satposes, zenith, dt, frequency), time[i])

        deltat = (jacobian)\-Rt[i]

        time[i+1]  = time[i] + deltat

    end
    
    Rt = Rt[1:iters]
    time = time[1:iters]
    time_measurment = time[end]
    
    return time_measurment

end

function tag_solve_TOA(all_sats_scaled, guess, zenith_angles, d, frequency)
    
    distance_scale = R_EARTH*(1e-3) # scales km to custom scale
    time_scale = 1/(C_LIGHT/R_EARTH) #scales seconds to custom scale
    c = 1 # speed of light

    #variables for the ionosphere time delay term
    Ip = zeros(4)
    OF = zeros(4)
    Iz = zeros(4)
    Ip_scaled = zeros(4)

    #Generate TEC Distribution for night TEC values. For 1 frequency
    mu = 8e16
    sigma = 3e16
    lb = 3e16
    ub = 13e16
    d_unscaled = Truncated(Normal(mu,sigma), lb, ub)



    n = 1000 # number of iterations
    #all_r0 = NaN*[zeros(4) for i = 1:n]

    #clock error distribution
    mu_t = 0
    sigma_t = 20e-9 #RMSE for Novatel GPS receiver
    lb_t = 1e-9
    ub_t = 30e-9
    d_t = Normal(mu_t, sigma_t)
    hi = 350e3 #in meters

    all_r0 = [zeros(4) for i = 1:n]
    sat1poses = zeros(4,n)
    sat2poses = zeros(4,n)
    sat3poses = zeros(4,n)
    sat4poses = zeros(4,n)

    all_sats_noise = zeros(4,4)

    #If there are 2 frequencies, add additional rows

    if frequency[2] != 0

        all_r0 = [zeros(8) for i = 1:n]
        sat5poses = zeros(4,n)
        sat6poses = zeros(4,n)
        sat7poses = zeros(4,n)
        sat8poses = zeros(4,n)

        all_sats_noise = zeros(8,4)

    end

    iters = 0

    #Monte Carlo Simulation
    for i in 1:n
        #Parameters for Armijo Line Search
        b = 0.01
        c= 0.5
        β = 10
        α = 1
        
        #create noise from random distribution
        #scaled to the variable (0.1m error for distance & 1e-11 for time)
        gpsnoise = randn(12)*(1e-4/distance_scale) #adding a 0.1 meter (10 cm) of noise 
        #clocknoise = randn(8)*(1e-11/time_scale) #original working (may not be accurate)
        #clocknoise = randn(8)*1e-9/time_scale # from Novatel resource ~20 ns accuracy
        clocknoise = rand(d_t,8)/time_scale

        TEC = rand(d_unscaled,4)

        #Calculating random TEC time delay
        for i=1:4
        
            OF[i] = (1-((R_EARTH*sind(zenith_angles[i]))/(R_EARTH+hi))^2)^(-1/2) #normal units (m and s)
            Ip[i] = ((40.3*TEC[i])/(frequency[1]^2))*OF[i] * 1e-3 #scale to km 
            Ip_scaled[i] = Ip[i]/distance_scale #scale to custom units
        
        end

        sat1_noise = [all_sats_scaled[1,1]+gpsnoise[1],all_sats_scaled[1,2] + gpsnoise[2],all_sats_scaled[1,3]+gpsnoise[3], all_sats_scaled[1,4] + clocknoise[1] + (Ip_scaled[1]/c)]
        sat2_noise =[all_sats_scaled[2,1]+gpsnoise[4],all_sats_scaled[2,2] + gpsnoise[5],all_sats_scaled[2,3]+gpsnoise[6], all_sats_scaled[2,4] + clocknoise[2] + (Ip_scaled[2]/c)]
        sat3_noise = [all_sats_scaled[3,1]+gpsnoise[7],all_sats_scaled[3,2] + gpsnoise[8],all_sats_scaled[3,3]+gpsnoise[9], all_sats_scaled[3,4] + clocknoise[3] + (Ip_scaled[3]/c)]
        sat4_noise = [all_sats_scaled[4,1]+gpsnoise[10],all_sats_scaled[4,2] + gpsnoise[11],all_sats_scaled[4,3]+gpsnoise[12], all_sats_scaled[4,4] + clocknoise[4] + (Ip_scaled[4]/c)]


        sat1poses[:,i] = sat1_noise
        sat2poses[:,i] = sat2_noise
        sat3poses[:,i] = sat3_noise
        sat4poses[:,i] = sat4_noise


        if frequency[2] != 0

            sat1_noise = [all_sats_scaled[1,1]+gpsnoise[1],all_sats_scaled[1,2] + gpsnoise[2],all_sats_scaled[1,3]+gpsnoise[3], all_sats_scaled[1,4] + clocknoise[1]]
            sat2_noise =[all_sats_scaled[2,1]+gpsnoise[4],all_sats_scaled[2,2] + gpsnoise[5],all_sats_scaled[2,3]+gpsnoise[6], all_sats_scaled[2,4] + clocknoise[2]]
            sat3_noise = [all_sats_scaled[3,1]+gpsnoise[7],all_sats_scaled[3,2] + gpsnoise[8],all_sats_scaled[3,3]+gpsnoise[9], all_sats_scaled[3,4] + clocknoise[3]]
            sat4_noise = [all_sats_scaled[4,1]+gpsnoise[10],all_sats_scaled[4,2] + gpsnoise[11],all_sats_scaled[4,3]+gpsnoise[12], all_sats_scaled[4,4] + clocknoise[4]]
            sat5_noise = [all_sats_scaled[1,1]+gpsnoise[1],all_sats_scaled[1,2] + gpsnoise[2],all_sats_scaled[1,3]+gpsnoise[3], all_sats_scaled[5,4] + clocknoise[5]]
            sat6_noise =[all_sats_scaled[2,1]+gpsnoise[4],all_sats_scaled[2,2] + gpsnoise[5],all_sats_scaled[2,3]+gpsnoise[6], all_sats_scaled[6,4] + clocknoise[6]]
            sat7_noise = [all_sats_scaled[3,1]+gpsnoise[7],all_sats_scaled[3,2] + gpsnoise[8],all_sats_scaled[3,3]+gpsnoise[9], all_sats_scaled[7,4] + clocknoise[7]]
            sat8_noise = [all_sats_scaled[4,1]+gpsnoise[10],all_sats_scaled[4,2] + gpsnoise[11],all_sats_scaled[4,3]+gpsnoise[12], all_sats_scaled[8,4] + clocknoise[8]]

            all_sats_noise[:,:] = vcat(sat1_noise',sat2_noise',sat3_noise',sat4_noise',sat5_noise',sat6_noise',sat7_noise',sat8_noise')
            
            #second frequency same 4 satellite config
            sat5poses[:,i] = sat5_noise
            sat6poses[:,i] = sat6_noise
            sat7poses[:,i] = sat7_noise
            sat8poses[:,i] = sat8_noise

        else

            all_sats_noise[:,:] = vcat(sat1_noise',sat2_noise',sat3_noise',sat4_noise')

        end
        
        X = NaN*[zeros(8) for i = 1:1000]
        R = NaN*[zeros(8) for i = 1:1000]

        if frequency[2] == 0
            #Centroid Guess
            rand_TEC = rand(d,4)
            X[1] = [guess[1]/distance_scale, guess[2]/distance_scale, guess[3]/distance_scale, 0.0001/time_scale, rand_TEC[1], rand_TEC[2], rand_TEC[3], rand_TEC[4]]

        else
            #Centroid guess
            X[1] = [guess[1]/distance_scale, guess[2]/distance_scale, guess[3]/distance_scale, 0.0001/time_scale, 3e-6, 3.5e-6, 2e-6, 4e-6]

        end

        #Initial residual

        R[1] = TOA_residual(X[1], all_sats_noise, zenith_angles, frequency)

        iters = 0
        
        #max iterations

        for k=1:1000
            
            R[k] = TOA_residual(X[k], all_sats_noise, zenith_angles, frequency)
            
            #println("this is residual: ", norm(R[k]))

            iters += 1
        
            #original tolerancec 1e-12
            
            if(norm(R[k]) < 1e-10)

                break

            end

            jacobian = ForwardDiff.jacobian(dx -> TOA_residual(dx, all_sats_noise, zenith_angles, frequency), X[k])

            conditionnum = cond(jacobian)
            
            #println("this is the jacobian: ", jacobian)
            
            #println("condition number: ", conditionnum)
            
            #println("This is R[k]: ", R[k])
            
            #println("this is TEC: ", TEC)
            #println("this is X[k]: ", X[k])
            
            #println("this is R[k]: ", norm(R[k]))
            
            #println("this is jacobian: ", jacobian)
            
            if frequency[2] == 0
                
                #only updating pose and clock offset
                deltax_1freq = (jacobian[1:4,1:4])\-R[k]

                deltax = vec([deltax_1freq zeros(4)])
            
                while norm(TOA_residual(X[k] + α*deltax, all_sats_noise, zenith_angles, frequency)) > norm(TOA_residual(X[k], all_sats_noise, zenith_angles, frequency) + b*α*jacobian[1:4,1:4]'*deltax[1:4])

                    α = c*α
                    #print("this is alpha: ", α)
                end

                X[k+1] = X[k] + α*deltax

            else
            
                deltax = (jacobian)\-R[k]
            
                while norm(TOA_residual(X[k] + α*deltax, all_sats_noise, zenith_angles, frequency)) > norm(TOA_residual(X[k], all_sats_noise, zenith_angles, frequency) + b*α*jacobian'*deltax)

                    α = c*α
                    #print("this is alpha: ", α)
                end

                X[k+1] = X[k] + α*deltax
            
            end
            
        end  
        
        R = R[1:iters]
        X = X[1:iters]

        all_r0[i] = X[end]

    end
    x_values = zeros(n)
    y_values = zeros(n)
    z_values = zeros(n)

    for j in 1:n

        x_values[j] = all_r0[j][1]
        y_values[j] = all_r0[j][2]
        z_values[j] = all_r0[j][3]

    end
    
    #mean calculation
    mean = sum(all_r0, dims = 1)/ n
    mean = mean[1]
    mean_rescaled = [mean[1:3]*distance_scale*1e3; mean[4]*time_scale; mean[5:8]] # in meters

    return mean_rescaled, all_r0, iters

end

function get_mean_TOA(combined_eci, i, d, r0_scaled, frequency, zenith_angles)

    distance_scale = R_EARTH*(1e-3) # scales km to custom scale
    time_scale = 1/(C_LIGHT/R_EARTH) #scales seconds to custom scale
    c = 1 # 1 distance_unit/ 1 time_unit

    TEC = rand(d,4) #vTEC for each satellite from custom distribution
    
    r0_TEC= vec([r0_scaled TEC]) #actual value

    sat_poses = [combined_eci[1,i] combined_eci[2,i] combined_eci[3,i];combined_eci[7,i] combined_eci[8,i] combined_eci[9,i]; combined_eci[13,i] combined_eci[14,i] combined_eci[15,i]; combined_eci[19,i] combined_eci[20,i] combined_eci[21,i]]*1e-3/distance_scale
     
    t = get_time(r0_TEC, sat_poses, zenith_angles, frequency)


    #In km then to custom units
    r1 = [combined_eci[1,i], combined_eci[2,i], combined_eci[3,i], 0]*1e-3/distance_scale
    r2 = [combined_eci[7,i], combined_eci[8,i], combined_eci[9,i], 0]*1e-3/distance_scale
    r3 = [combined_eci[13,i],combined_eci[14,i], combined_eci[15,i], 0]*1e-3/distance_scale
    r4 = [combined_eci[19,i], combined_eci[20,i], combined_eci[21,i], 0]*1e-3/distance_scale
    
    r1[4] = t[1] 
    r2[4] = t[2]
    r3[4] = t[3]
    r4[4] = t[4]

    all_sats_scaled = vcat(r1',r2',r3',r4')
    
    if frequency[2] != 0
        
        #measurments at 2nd frequency
        r5 = [combined_eci[1,i], combined_eci[2,i], combined_eci[3,i], 0]*1e-3/distance_scale
        r6 = [combined_eci[7,i], combined_eci[8,i], combined_eci[9,i], 0]*1e-3/distance_scale
        r7 = [combined_eci[13,i],combined_eci[14,i], combined_eci[15,i], 0]*1e-3/distance_scale
        r8 = [combined_eci[19,i], combined_eci[20,i], combined_eci[21,i], 0]*1e-3/distance_scale
        r5[4] = t[5] 
        r6[4] = t[6]
        r7[4] = t[7]
        r8[4] = t[8]

        all_sats_scaled = vcat(r1',r2',r3',r4', r5', r6', r7', r8')
    end
    
    #Generate a guess at the centroid of all 4 sats on the surface of the Earth. Scaled to km
    centroid_guess = [(combined_eci[1,i]+combined_eci[7,i]+combined_eci[13,i]+combined_eci[19,i])/4, (combined_eci[2,i]+combined_eci[8,i]+combined_eci[14,i]+combined_eci[20,i])/4, (combined_eci[3,i]+combined_eci[9,i]+combined_eci[15,i]+combined_eci[21,i])/4] 
    onearth = sECEFtoGEOC(centroid_guess, use_degrees = true)
    geodetic = [onearth[1], onearth[2], 0]

    xyz_guess = sGEOCtoECEF(geodetic, use_degrees = true)*1e-3 #change to km
    
    mean_rescaled, all_r0, iters = tag_solve_TOA(all_sats_scaled, xyz_guess, zenith_angles, d, frequency)

    return mean_rescaled, all_r0, iters
end

function find_covariance(all_r0, mean_rescaled)

    distance_scale = R_EARTH*(1e-3) # scales km to custom scale
    time_scale = 1/(C_LIGHT/R_EARTH) #scales seconds to custom scale
    c = 1 # 1 distance_unit/ 1 time_unit

    n=1000

    all_r0_scaled = zeros(n,8)

    #Rescale back to units (m and s)
    for i in 1:n
    
        all_r0_scaled[i,:] = [all_r0[i][1:3]*distance_scale*1e3; all_r0[i][4]*time_scale; all_r0[i][5:8]]
        
    end
    
    totalsum = zeros(8,8)

    for i in 1:n
    
        value = all_r0_scaled[i,:] - mean_rescaled
        matrixvalue = value*transpose(value)
        totalsum += matrixvalue
    
    end

    covariancematrix = totalsum/n
    
    PDOP = sqrt(covariancematrix[1,1] + covariancematrix[2,2] + covariancematrix[3,3])
    TDOP = sqrt(covariancematrix[4,4])
    #TEC_cov = sqrt(covariancematrix[5,5] + covariancematrix[6,6] + covariancematrix[7,7] + covariancematrix[8,8])
    
    return PDOP

end

export test

end # module
