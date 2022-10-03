import SatelliteNavigation as SN

#Doppler Functions

function doppler_residual(x, sat_pose_f, time, z, zdot, frequency)

    #Doppler Custom Units
    distance_scale = R_EARTH*(1e-3) # scales km to custom scale working

    time_scale = (1/(C_LIGHT/R_EARTH))*1e4 #scales seconds to custom scale

    velocity_scale = distance_scale/time_scale #scales velocities to custom scales

    frequency_scale = 1

    c = 1*1e4 #km/s

    h = 350e3 #ionosphere thin shell approx height in meters
    
    omega = OMEGA_EARTH*time_scale #transform to custom scale
        
    omegahat = [0 -omega 0; omega 0 0; 0 0 0]
    
    A1 = [cos(omega*time[1]) sin(omega*time[1]) 0;
          -sin(omega*time[1]) cos(omega*time[1]) 0;
          0 0 1];
    A2 = [cos(omega*time[2]) sin(omega*time[2]) 0;
          -sin(omega*time[2]) cos(omega*time[2]) 0;
          0 0 1];
    A3 = [cos(omega*time[3]) sin(omega*time[3]) 0;
          -sin(omega*time[3]) cos(omega*time[3]) 0;
          0 0 1];
    A4 = [cos(omega*time[4]) sin(omega*time[4]) 0;
          -sin(omega*time[4]) cos(omega*time[4]) 0;
          0 0 1];

    #if we have 2 frequencies, we are estimating out the TEC values as part of the state
    if frequency[2] != 0

        #bdot is a velocity. therefore it is scaled
        #Ionosphere Errors at frequency 1
        Idot1 = (((40.3*x[5])*((R_EARTH*sind(z[1]))*(R_EARTH*cosd(z[1])*zdot[1])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[1]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot2 = (((40.3*x[6])*((R_EARTH*sind(z[2]))*(R_EARTH*cosd(z[2])*zdot[2])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[2]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot3 = (((40.3*x[7])*((R_EARTH*sind(z[3]))*(R_EARTH*cosd(z[3])*zdot[3])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[3]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot4 = (((40.3*x[8])*((R_EARTH*sind(z[4]))*(R_EARTH*cosd(z[4])*zdot[4])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[4]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale

        #Ionosphere Errors at frequency 2
        Idot5 = (((40.3*x[5])*((R_EARTH*sind(z[1]))*(R_EARTH*cosd(z[1])*zdot[1])/(R_EARTH + h)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[1]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot6 = (((40.3*x[6])*((R_EARTH*sind(z[2]))*(R_EARTH*cosd(z[2])*zdot[2])/(R_EARTH + h)^2))/((frequency[2])*(1-((R_EARTH*sind(z[2]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot7 = (((40.3*x[7])*((R_EARTH*sind(z[3]))*(R_EARTH*cosd(z[3])*zdot[3])/(R_EARTH + h)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[3]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot8 = (((40.3*x[8])*((R_EARTH*sind(z[4]))*(R_EARTH*cosd(z[4])*zdot[4])/(R_EARTH + h)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[4]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale


    else
        #ionosphere effect is added later for the 1 frequency case
        Idot1 = 0
        Idot2 = 0
        Idot3 = 0
        Idot4 = 0
    end

    #residuals to estimate a frequency offset bdot
    res1 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A1*sat_pose_f[1,1:3]))*(-2*x[1:3]'*A1*omegahat*sat_pose_f[1,1:3] - 2*x[1:3]'*A1*sat_pose_f[1,4:6]+sat_pose_f[1,1:3]'*sat_pose_f[1,4:6]+sat_pose_f[1,4:6]'*sat_pose_f[1,1:3])+Idot1 + x[4]) - sat_pose_f[1,7]
    res2 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A2*sat_pose_f[2,1:3]))*(-2*x[1:3]'*A2*omegahat*sat_pose_f[2,1:3] - 2*x[1:3]'*A2*sat_pose_f[2,4:6]+sat_pose_f[2,1:3]'*sat_pose_f[2,4:6]+sat_pose_f[2,4:6]'*sat_pose_f[2,1:3]) +Idot2 + x[4]) - sat_pose_f[2,7]
    res3 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A3*sat_pose_f[3,1:3]))*(-2*x[1:3]'*A3*omegahat*sat_pose_f[3,1:3] - 2*x[1:3]'*A3*sat_pose_f[3,4:6]+sat_pose_f[3,1:3]'*sat_pose_f[3,4:6]+sat_pose_f[3,4:6]'*sat_pose_f[3,1:3]) +Idot3 +  x[4]) - sat_pose_f[3,7]
    res4 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A4*sat_pose_f[4,1:3]))*(-2*x[1:3]'*A4*omegahat*sat_pose_f[4,1:3] - 2*x[1:3]'*A4*sat_pose_f[4,4:6]+sat_pose_f[4,1:3]'*sat_pose_f[4,4:6]+sat_pose_f[4,4:6]'*sat_pose_f[4,1:3]) + Idot4 + x[4]) - sat_pose_f[4,7]
    
    residuals = [res1, res2, res3, res4]

    if frequency[2] != 0
        #Measurments at frequency 2
        res5 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A1*sat_pose_f[1,1:3]))*(-2*x[1:3]'*A1*omegahat*sat_pose_f[1,1:3] - 2*x[1:3]'*A1*sat_pose_f[1,4:6]+sat_pose_f[1,1:3]'*sat_pose_f[1,4:6]+sat_pose_f[1,4:6]'*sat_pose_f[1,1:3]) + Idot5 + x[4]) - sat_pose_f[5,7]
        res6 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A2*sat_pose_f[2,1:3]))*(-2*x[1:3]'*A2*omegahat*sat_pose_f[2,1:3] - 2*x[1:3]'*A2*sat_pose_f[2,4:6]+sat_pose_f[2,1:3]'*sat_pose_f[2,4:6]+sat_pose_f[2,4:6]'*sat_pose_f[2,1:3]) + Idot6 + x[4]) - sat_pose_f[6,7]
        res7 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A3*sat_pose_f[3,1:3]))*(-2*x[1:3]'*A3*omegahat*sat_pose_f[3,1:3] - 2*x[1:3]'*A3*sat_pose_f[3,4:6]+sat_pose_f[3,1:3]'*sat_pose_f[3,4:6]+sat_pose_f[3,4:6]'*sat_pose_f[3,1:3]) + Idot7 + x[4]) - sat_pose_f[7,7]
        res8 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(x[1:3] - A4*sat_pose_f[4,1:3]))*(-2*x[1:3]'*A4*omegahat*sat_pose_f[4,1:3] - 2*x[1:3]'*A4*sat_pose_f[4,4:6]+sat_pose_f[4,1:3]'*sat_pose_f[4,4:6]+sat_pose_f[4,4:6]'*sat_pose_f[4,1:3]) + Idot8 + x[4]) - sat_pose_f[8,7]
        
        residuals = [res1, res2, res3, res4, res5, res6, res7, res8]

    end

    return residuals
   
end

#pose1-4 are satellite locations and velocities
#r0 = [xyz of tag, emmitted frequency]

function doppler_measurment(r0, sat_poses, time, z, zdot, frequency)


    #Doppler Custom Units
    distance_scale = R_EARTH*(1e-3) # scales km to custom scale working

    time_scale = (1/(C_LIGHT/R_EARTH))*1e4 #scales seconds to custom scale

    velocity_scale = distance_scale/time_scale #scales velocities to custom scales

    frequency_scale = 1

    c = 1*1e4 #km/s

    h = 350e3 #ionosphere thin shell approx height in meters

    #Make sure to scale all positions
    
    omega = OMEGA_EARTH*time_scale #change to custom time scale
    
    #wavelength = 0.7494*1e-3/distance_scale #for a nominal frequency of 400 MHz
    
    omegahat = [0 -omega 0; omega 0 0; 0 0 0]
    
    A1 = [cos(omega*time[1]) sin(omega*time[1]) 0;
          -sin(omega*time[1]) cos(omega*time[1]) 0;
          0 0 1];
    A2 = [cos(omega*time[2]) sin(omega*time[2]) 0;
          -sin(omega*time[2]) cos(omega*time[2]) 0;
          0 0 1];
    A3 = [cos(omega*time[3]) sin(omega*time[3]) 0;
          -sin(omega*time[3]) cos(omega*time[3]) 0;
          0 0 1];
    A4 = [cos(omega*time[4]) sin(omega*time[4]) 0;
          -sin(omega*time[4]) cos(omega*time[4]) 0;
          0 0 1];

    if frequency[2] != 0

    #Ionosphere Errors at frequency 1
        Idot1 = (((40.3*r0[5])*((R_EARTH*sind(z[1]))*(R_EARTH*cosd(z[1])*zdot[1])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[1]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot2 = (((40.3*r0[6])*((R_EARTH*sind(z[2]))*(R_EARTH*cosd(z[2])*zdot[2])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[2]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot3 = (((40.3*r0[7])*((R_EARTH*sind(z[3]))*(R_EARTH*cosd(z[3])*zdot[3])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[3]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot4 = (((40.3*r0[8])*((R_EARTH*sind(z[4]))*(R_EARTH*cosd(z[4])*zdot[4])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[4]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale

        #Ionosphere Errors at frequency 1
        Idot5 = (((40.3*r0[5])*((R_EARTH*sind(z[1]))*(R_EARTH*cosd(z[1])*zdot[1])/(R_EARTH + h)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[1]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot6 = (((40.3*r0[6])*((R_EARTH*sind(z[2]))*(R_EARTH*cosd(z[2])*zdot[2])/(R_EARTH + h)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[2]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot7 = (((40.3*r0[7])*((R_EARTH*sind(z[3]))*(R_EARTH*cosd(z[3])*zdot[3])/(R_EARTH + h)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[3]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot8 = (((40.3*r0[8])*((R_EARTH*sind(z[4]))*(R_EARTH*cosd(z[4])*zdot[4])/(R_EARTH + h)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[4]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
    else
        Idot1 = 0
        Idot2 = 0
        Idot3 = 0
        Idot4 = 0

    end
    #measurments for estimating a transmitting frequency
    #deltaf1 = (r0[4]/c)*0.5*(1/norm(r0[1:3] - A1*pose1[1:3]))*(-2*r0[1:3]'*A1*omegahat*pose1[1:3] - 2*r0[1:3]'*A1*pose1[4:6]+pose1[1:3]'*pose1[4:6]+pose1[4:6]'*pose1[1:3])
    #deltaf2 = (r0[4]/c)*0.5*(1/norm(r0[1:3] - A2*pose2[1:3]))*(-2*r0[1:3]'*A2*omegahat*pose2[1:3] - 2*r0[1:3]'*A2*pose2[4:6]+pose2[1:3]'*pose2[4:6]+pose2[4:6]'*pose2[1:3])
    #deltaf3 = (r0[4]/c)*0.5*(1/norm(r0[1:3] - A3*pose3[1:3]))*(-2*r0[1:3]'*A3*omegahat*pose3[1:3] - 2*r0[1:3]'*A3*pose3[4:6]+pose3[1:3]'*pose3[4:6]+pose3[4:6]'*pose3[1:3])
    #deltaf4 = (r0[4]/c)*0.5*(1/norm(r0[1:3] - A4*pose4[1:3]))*(-2*r0[1:3]'*A4*omegahat*pose4[1:3] - 2*r0[1:3]'*A4*pose4[4:6]+pose4[1:3]'*pose4[4:6]+pose4[4:6]'*pose4[1:3])

    #measurments to estimate out a bdot value (frequency offset)
    deltaf1 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A1*sat_poses[1, 1:3]))*(-2*r0[1:3]'*A1*omegahat*sat_poses[1, 1:3] - 2*r0[1:3]'*A1*sat_poses[1, 4:6]+sat_poses[1, 1:3]'*sat_poses[1, 4:6]+sat_poses[1, 4:6]'*sat_poses[1, 1:3]) + Idot1 + r0[4])
    deltaf2 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A2*sat_poses[2, 1:3]))*(-2*r0[1:3]'*A2*omegahat*sat_poses[2, 1:3] - 2*r0[1:3]'*A2*sat_poses[2, 4:6]+sat_poses[2, 1:3]'*sat_poses[2, 4:6]+sat_poses[2, 4:6]'*sat_poses[2, 1:3]) + Idot2 + r0[4])
    deltaf3 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A3*sat_poses[3, 1:3]))*(-2*r0[1:3]'*A3*omegahat*sat_poses[3, 1:3] - 2*r0[1:3]'*A3*sat_poses[3, 4:6]+sat_poses[3, 1:3]'*sat_poses[3, 4:6]+sat_poses[3, 4:6]'*sat_poses[3, 1:3]) + Idot3 + r0[4])
    deltaf4 = (frequency[1]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A4*sat_poses[4, 1:3]))*(-2*r0[1:3]'*A4*omegahat*sat_poses[4, 1:3] - 2*r0[1:3]'*A4*sat_poses[4, 4:6]+sat_poses[4, 1:3]'*sat_poses[4, 4:6]+sat_poses[4, 4:6]'*sat_poses[4, 1:3]) + Idot4 + r0[4])
    
    deltaf = [deltaf1, deltaf2, deltaf3, deltaf4]

    if frequency[2] != 0 

    #Measurments at frequency 2
    deltaf5 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A1*sat_poses[1, 1:3]))*(-2*r0[1:3]'*A1*omegahat*sat_poses[1, 1:3] - 2*r0[1:3]'*A1*sat_poses[1, 4:6]+sat_poses[1, 1:3]'*sat_poses[1, 4:6]+sat_poses[1, 4:6]'*sat_poses[1, 1:3]) + Idot5 + r0[4])
    deltaf6 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A2*sat_poses[2, 1:3]))*(-2*r0[1:3]'*A2*omegahat*sat_poses[2, 1:3] - 2*r0[1:3]'*A2*sat_poses[2, 4:6]+sat_poses[2, 1:3]'*sat_poses[2, 4:6]+sat_poses[2, 4:6]'*sat_poses[2, 1:3]) + Idot6 + r0[4])
    deltaf7 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A3*sat_poses[3, 1:3]))*(-2*r0[1:3]'*A3*omegahat*sat_poses[3, 1:3] - 2*r0[1:3]'*A3*sat_poses[3, 4:6]+sat_poses[3, 1:3]'*sat_poses[3, 4:6]+sat_poses[3, 4:6]'*sat_poses[3, 1:3]) + Idot7 + r0[4])
    deltaf8 = (frequency[2]*frequency_scale/c)*(0.5*(1/norm(r0[1:3] - A4*sat_poses[4, 1:3]))*(-2*r0[1:3]'*A4*omegahat*sat_poses[4, 1:3] - 2*r0[1:3]'*A4*sat_poses[4, 4:6]+sat_poses[4, 1:3]'*sat_poses[4, 4:6]+sat_poses[4, 4:6]'*sat_poses[4, 1:3]) + Idot8 + r0[4])
    
    deltaf = [deltaf1, deltaf2, deltaf3, deltaf4, deltaf5, deltaf6, deltaf7, deltaf8]
    
    end

    return deltaf
end


function tag_solve_Doppler(all_sats_scaled, guess, zenith_angles, time, zdot, frequency, d, r0_scaled) # remember to scale the guess to custom units

    #Doppler Custom Units
    distance_scale = R_EARTH*(1e-3) # scales km to custom scale working

    time_scale = (1/(C_LIGHT/R_EARTH))*1e4 #scales seconds to custom scale

    velocity_scale = distance_scale/time_scale #scales velocities to custom scales

    frequency_scale = 1

    c = 1*1e4 #km/s

    h = 350e3 #ionosphere thin shell approx height in meters

    #works well for double frequency
    #diff = 50e3*1e-3/distance_scale 

    diff = 25e3*1e-3/distance_scale 


    n = 1000 # number of iterations
    #all_r0 = NaN*[zeros(4) for i = 1:n]
    if frequency[2] != 0

        all_r0 = [zeros(8) for i = 1:n] #obtain all the tag positions and frequency offset

        sat1poses = zeros(7,n)
        sat2poses = zeros(7,n)
        sat3poses = zeros(7,n)
        sat4poses = zeros(7,n)
        sat5poses = zeros(7,n)
        sat6poses = zeros(7,n)
        sat7poses = zeros(7,n)
        sat8poses = zeros(7,n)

        all_sats_noise = zeros(8,4)
    
    else

        all_r0 = [zeros(4) for i = 1:n] #obtain all the tag positions
        
        sat1poses = zeros(7,n)
        sat2poses = zeros(7,n)
        sat3poses = zeros(7,n)
        sat4poses = zeros(7,n)

        all_sats_noise = zeros(4,4)

    end

    iters = 0

     #Monte Carlo Simulation
    for i in 1:n
        #Parameters for line search
        #b = 0.01
        b = 0.1
        c_=0.5
        #β = 10
        β = 1.0
        α = 1
        
        #create noise from normal distribution
        
        #scaled to the variable (0.1m error for distance & 1e-11 for time)
        gpsnoise = randn(12)*(0.1*1e-3/distance_scale) #0.1 m switch to km then to custom scale
        velocitynoise = randn(12)*(0.01*1e-3/velocity_scale)# 0.1 m/s
        doppler_noise = randn(12)* 1/frequency_scale #0.001 ~ 500's of error
        
        #doppler_noise = randn(12)* 0/frequency_scale #0.001 ~ 500's of error
        
        TEC = rand(d,4) #vTEC for each satellite from custom distribution
        

        #Calculating random TEC time delay
        #for i=1:4
        
        #    OF[i] = (1-((R_EARTH*sind(zenith_angles[i]))/(R_EARTH+hi))^2)^(-1/2) #normal units (m and s)
        #    Ip[i] = ((40.3*TEC[i])/(f^2))*OF[i] * 1e-3 #scale to km to use the custom unit
        #    Ip_scaled[i] = Ip[i]/distance_scale #scale to custom units
        
        #end


        Idot1 = (((40.3*TEC[1])*((R_EARTH*sind(zenith_angles[1]))*(R_EARTH*cosd(zenith_angles[1])*zdot[1])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[1]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot2 = (((40.3*TEC[2])*((R_EARTH*sind(zenith_angles[2]))*(R_EARTH*cosd(zenith_angles[2])*zdot[2])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[2]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot3 = (((40.3*TEC[3])*((R_EARTH*sind(zenith_angles[3]))*(R_EARTH*cosd(zenith_angles[3])*zdot[3])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[3]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        Idot4 = (((40.3*TEC[4])*((R_EARTH*sind(zenith_angles[4]))*(R_EARTH*cosd(zenith_angles[4])*zdot[4])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[4]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale

        #print("these are zenith angles: ", zenith_angles)
        #Idot1 = (((40.3*TEC[1])*((R_EARTH*sind(zenith_angles[1]))*(R_EARTH*cosd(zenith_angles[1])*zdot[1])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[1]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        #Idot2 = (((40.3*TEC[2])*((R_EARTH*sind(zenith_angles[2]))*(R_EARTH*cosd(zenith_angles[2])*zdot[2])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[2]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        #Idot3 = (((40.3*TEC[3])*((R_EARTH*sind(zenith_angles[3]))*(R_EARTH*cosd(zenith_angles[3])*zdot[3])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[3]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        #Idot4 = (((40.3*TEC[4])*((R_EARTH*sind(zenith_angles[4]))*(R_EARTH*cosd(zenith_angles[4])*zdot[4])/(R_EARTH + h)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(zenith_angles[4]))/(R_EARTH+h))^2)^1.5))*1e20*1e-3/velocity_scale
        
        iono1 = (frequency[1]/c)*Idot1
        iono2 = (frequency[1]/c)*Idot2
        iono3 = (frequency[1]/c)*Idot3
        iono4 = (frequency[1]/c)*Idot4
        
        #print("this is f: ", frequency[1])
        #print("this is c: ", c)
        #print("this is Idot1: ", Idot1)

        scaled_iono1 = iono1/distance_scale*velocity_scale
        
        #println("this is scaled_iono1: ", scaled_iono1) #in Hz
        
        sat1_noise = [all_sats_scaled[1,1]+gpsnoise[1],all_sats_scaled[1,2] + gpsnoise[2],all_sats_scaled[1,3]+gpsnoise[3], all_sats_scaled[1,4]+velocitynoise[1], all_sats_scaled[1,5]+velocitynoise[2],all_sats_scaled[1,6]+velocitynoise[3], all_sats_scaled[1,7]+iono1+doppler_noise[1]]
        sat2_noise =[all_sats_scaled[2,1]+gpsnoise[4],all_sats_scaled[2,2] + gpsnoise[5],all_sats_scaled[2,3]+gpsnoise[6], all_sats_scaled[2,4]+velocitynoise[4], all_sats_scaled[2,5]+velocitynoise[5],all_sats_scaled[2,6]+velocitynoise[6], all_sats_scaled[2,7]+iono2+doppler_noise[2]]
        sat3_noise = [all_sats_scaled[3,1]+gpsnoise[7],all_sats_scaled[3,2] + gpsnoise[8],all_sats_scaled[3,3]+gpsnoise[9], all_sats_scaled[3,4]+velocitynoise[7], all_sats_scaled[3,5]+velocitynoise[8],all_sats_scaled[3,6]+velocitynoise[9], all_sats_scaled[3,7]+iono3+doppler_noise[3]]
        sat4_noise = [all_sats_scaled[4,1]+gpsnoise[10],all_sats_scaled[4,2] + gpsnoise[11],all_sats_scaled[4,3]+gpsnoise[12], all_sats_scaled[4,4]+velocitynoise[10], all_sats_scaled[4,5]+velocitynoise[11],all_sats_scaled[4,6]+velocitynoise[12], all_sats_scaled[4,7]+iono4+doppler_noise[4]]
        
        sat1poses[:,i] = sat1_noise
        sat2poses[:,i] = sat2_noise
        sat3poses[:,i] = sat3_noise
        sat4poses[:,i] = sat4_noise

        if frequency[2] != 0

            sat1_noise = [all_sats_scaled[1,1]+gpsnoise[1],all_sats_scaled[1,2] + gpsnoise[2],all_sats_scaled[1,3]+gpsnoise[3], all_sats_scaled[1,4]+velocitynoise[1], all_sats_scaled[1,5]+velocitynoise[2],all_sats_scaled[1,6]+velocitynoise[3], all_sats_scaled[1,7] + doppler_noise[1]]
            sat2_noise =[all_sats_scaled[2,1]+gpsnoise[4],all_sats_scaled[2,2] + gpsnoise[5],all_sats_scaled[2,3]+gpsnoise[6], all_sats_scaled[2,4]+velocitynoise[4], all_sats_scaled[2,5]+velocitynoise[5],all_sats_scaled[2,6]+velocitynoise[6], all_sats_scaled[2,7] + doppler_noise[2]]
            sat3_noise = [all_sats_scaled[3,1]+gpsnoise[7],all_sats_scaled[3,2] + gpsnoise[8],all_sats_scaled[3,3]+gpsnoise[9], all_sats_scaled[3,4]+velocitynoise[7], all_sats_scaled[3,5]+velocitynoise[8],all_sats_scaled[3,6]+velocitynoise[9], all_sats_scaled[3,7] + doppler_noise[3]]
            sat4_noise = [all_sats_scaled[4,1]+gpsnoise[10],all_sats_scaled[4,2] + gpsnoise[11],all_sats_scaled[4,3]+gpsnoise[12], all_sats_scaled[4,4]+velocitynoise[10], all_sats_scaled[4,5]+velocitynoise[11],all_sats_scaled[4,6]+velocitynoise[12], all_sats_scaled[4,7] + doppler_noise[4]]
            
            sat5_noise = [all_sats_scaled[1,1]+gpsnoise[1],all_sats_scaled[1,2] + gpsnoise[2],all_sats_scaled[1,3]+gpsnoise[3], all_sats_scaled[1,4]+velocitynoise[1], all_sats_scaled[1,5]+velocitynoise[2],all_sats_scaled[1,6]+velocitynoise[3], all_sats_scaled[5,7] + doppler_noise[5]]
            sat6_noise =[all_sats_scaled[2,1]+gpsnoise[4],all_sats_scaled[2,2] + gpsnoise[5],all_sats_scaled[2,3]+gpsnoise[6], all_sats_scaled[2,4]+velocitynoise[4], all_sats_scaled[2,5]+velocitynoise[5],all_sats_scaled[2,6]+velocitynoise[6], all_sats_scaled[6,7] + doppler_noise[6]]
            sat7_noise = [all_sats_scaled[3,1]+gpsnoise[7],all_sats_scaled[3,2] + gpsnoise[8],all_sats_scaled[3,3]+gpsnoise[9], all_sats_scaled[3,4]+velocitynoise[7], all_sats_scaled[3,5]+velocitynoise[8],all_sats_scaled[3,6]+velocitynoise[9], all_sats_scaled[7,7] + doppler_noise[7]]
            sat8_noise = [all_sats_scaled[4,1]+gpsnoise[10],all_sats_scaled[4,2] + gpsnoise[11],all_sats_scaled[4,3]+gpsnoise[12], all_sats_scaled[4,4]+velocitynoise[10], all_sats_scaled[4,5]+velocitynoise[11],all_sats_scaled[4,6]+velocitynoise[12], all_sats_scaled[8,7] + doppler_noise[8]]

            sat1poses[:,i] = sat1_noise
            sat2poses[:,i] = sat2_noise
            sat3poses[:,i] = sat3_noise
            sat4poses[:,i] = sat4_noise
            sat5poses[:,i] = sat5_noise
            sat6poses[:,i] = sat6_noise
            sat7poses[:,i] = sat7_noise
            sat8poses[:,i] = sat8_noise

            all_sats_noise = vcat(sat1_noise',sat2_noise',sat3_noise',sat4_noise', sat5_noise',sat6_noise',sat7_noise',sat8_noise')

        else

            all_sats_noise = vcat(sat1_noise',sat2_noise',sat3_noise',sat4_noise')

        end

        
        X = NaN*[zeros(8) for i = 1:n]
        R = NaN*[zeros(8) for i = 1:n]
        
        #Centroid guess. Working for estimating transmit frequency
        #X[1] = [guess[1]/distance_scale, guess[2]/distance_scale, guess[3]/distance_scale, 410e6/frequency_scale]
        
        
        #working
        X[1] = [r0_scaled[1]+diff, r0_scaled[2]+diff, r0_scaled[3]+diff, 1.002*1e4*1e-3/velocity_scale, 3.1e-4, 4.2e-4, 3.2e-4, 4.7e-4]

        if frequency[2] == 0
            X[1] = [r0_scaled[1]+diff, r0_scaled[2]+diff, r0_scaled[3]+diff, 1.002*1e4*1e-3/velocity_scale, 0,0,0,0]
        end
        #New guess to estimate bdot
        #X[1] = [guess[1]/distance_scale, guess[2]/distance_scale, guess[3]/distance_scale, 1.002*1e4*1e-3/velocity_scale]

        #X[1] = r0_scaled
        
        #println("this is initial guess: ", X[1])
        #Initial residual
        #R[1] = doppler_residual(X[1], all_sats_noise, time)

        R[1] = doppler_residual(X[1], all_sats_noise, time, zenith_angles, zdot, frequency)

        iters = 0

        for k=1:1000
            #R[k] = doppler_residual(X[k], all_sats_noise, time) #calculate residual

            R[k] = doppler_residual(X[k], all_sats_noise, time, zenith_angles, zdot, frequency)

            iters += 1

            #if(norm(R[k]) < 1e-6)
            
            if(norm(R[k]) < 1e-6)
                
                break

            end

            #jacobian = ForwardDiff.jacobian(dx -> doppler_residual(dx, all_sats_noise, time), X[k])

            jacobian = ForwardDiff.jacobian(dx -> doppler_residual(dx, all_sats_noise, time, zenith_angles, zdot, frequency), X[k])


            if frequency[2] != 0

                if(rank(jacobian) != 8) #if not full rank, move on to next iteration
                    
                    break

                end

            else

                if(rank(jacobian) != 4) #if not full rank, move on to next iteration
                    
                    break

                end
            

            end

            conditionnum = cond(jacobian)
            
            #println("this is residual: ", norm(R[k]))
            
            #println("this is the condition number: ", conditionnum)
            
            #println("this is the jacobian: ", jacobian)
            
            #println("this is X[k]: ", X[k])
            
            if iters > 50
                break
            end
            

            if frequency[2] == 0

                #only updating pose and frequency offset
                deltax_1freq = (jacobian[1:4,1:4])\-R[k]

                deltax = vec([deltax_1freq zeros(4)])
            
                #while norm(TOA_residual(X[k] + α*deltax, all_sats_noise, zenith_angles, frequency)) > norm(TOA_residual(X[k], all_sats_noise, zenith_angles, frequency) + b*α*jacobian[1:4,1:4]'*deltax[1:4])
                while norm(doppler_residual(X[k] + α*deltax, all_sats_noise, time, zenith_angles, zdot, frequency)) > norm(doppler_residual(X[k] + α*deltax, all_sats_noise, time, zenith_angles, zdot, frequency) + b*α*jacobian[1:4,1:4]'*deltax[1:4])
                    α = c*α
                    #print("this is alpha: ", α)
                end

                X[k+1] = X[k] + α*deltax

            else

                #2 frequencies use the whole state

                deltax = (jacobian)\-R[k]

                #while norm(doppler_residual(X[k] + α*deltax, all_sats_noise, time)) > norm(doppler_residual(X[k], all_sats_noise, time) + b*α*jacobian'*deltax)
                while norm(doppler_residual(X[k] + α*deltax, all_sats_noise, time, zenith_angles, zdot, frequency)) > norm(doppler_residual(X[k] + α*deltax, all_sats_noise, time, zenith_angles, zdot, frequency) + b*α*jacobian'*deltax)

                    α = c_*α

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
    #for estimating frequency
    #mean_rescaled = [mean[1:3]*distance_scale*1e3; mean[4]*frequency_scale] # in meters
    
    #For estimating bdot
    mean_rescaled = [mean[1:3]*distance_scale*1e3; mean[4]*1e3*velocity_scale; mean[5:8]] # in meters

    return mean_rescaled, all_r0, iters
end


function get_mean_Doppler(combined_eci, i, d, r0_scaled, frequency, zenith_angles, satnum, tag)


    #Doppler Custom Units
    distance_scale = R_EARTH*(1e-3) # scales km to custom scale working

    time_scale = (1/(C_LIGHT/R_EARTH))*1e4 #scales seconds to custom scale

    velocity_scale = distance_scale/time_scale #scales velocities to custom scales

    frequency_scale = 1

    c = 1*1e4 #km/s


    TEC = rand(d,4) #vTEC for each satellite from custom distribution
    
    r0_TEC= vec([r0_scaled TEC]) #actual value

    #print("sent zenith angles: ", zenith_angles)
    #calculate zdot 
    satpositions1 = [combined_eci[1,i] combined_eci[2,i] combined_eci[3,i];combined_eci[7,i] combined_eci[8,i] combined_eci[9,i]; combined_eci[13,i] combined_eci[14,i] combined_eci[15,i]; combined_eci[19,i] combined_eci[20,i] combined_eci[21,i]]
    satpositions2 = [combined_eci[1,i+1] combined_eci[2,i+1] combined_eci[3,i+1];combined_eci[7,i+1] combined_eci[8,i+1] combined_eci[9,i+1]; combined_eci[13,i+1] combined_eci[14,i+1] combined_eci[15,i+1]; combined_eci[19,i+1] combined_eci[20,i+1] combined_eci[21,i+1]]

    zenith_angles = SN.zenith(satpositions1, tag, satnum)
    zenith_angles2 = SN.zenith(satpositions2, tag, satnum)
 
    #print("newly calculated zenith angles: ", zenith_angles)
    zdot = zenith_angles - zenith_angles2


    r1 = [combined_eci[1,i], combined_eci[2,i], combined_eci[3,i]]*1e-3/distance_scale
    r2 = [combined_eci[7,i], combined_eci[8,i], combined_eci[9,i]]*1e-3/distance_scale
    r3 = [combined_eci[13,i], combined_eci[14,i], combined_eci[15,i]]*1e-3/distance_scale
    r4 = [combined_eci[19,i], combined_eci[20,i], combined_eci[21,i]]*1e-3/distance_scale

    v1 = [combined_eci[4,i], combined_eci[5,i], combined_eci[6,i]]*1e-3/velocity_scale
    v2 = [combined_eci[10,i], combined_eci[11,i], combined_eci[12,i]]*1e-3/velocity_scale
    v3 = [combined_eci[16,i], combined_eci[17,i], combined_eci[18,i]]*1e-3/velocity_scale
    v4 = [combined_eci[22,i], combined_eci[23,i], combined_eci[24,i]]*1e-3/velocity_scale

    pose1 = [r1' v1']
    pose2 = [r2' v2']
    pose3 = [r3' v3']
    pose4 = [r4' v4']

    sat_poses = [pose1; pose2; pose3; pose4]

    #sat_poses = [combined_eci[1:6,i]'; combined_eci[7:12, i]'; combined_eci[13:18, i]'; combined_eci[19:24, i]']
    #sat_poses = [combined_eci[1,i] combined_eci[2,i] combined_eci[3,i];combined_eci[7,i] combined_eci[8,i] combined_eci[9,i]; combined_eci[13,i] combined_eci[14,i] combined_eci[15,i]; combined_eci[19,i] combined_eci[20,i] combined_eci[21,i]]*1e-3/distance_scale
    
    #arbitrary term
    time = [0.006, 0.006, 0.006, 0.008]/time_scale # assume a fixed time for the rotation matrix


    #implement sat_poses into doppler_measurment
    deltaf = doppler_measurment(r0_TEC, sat_poses , time, zenith_angles, zdot, frequency)


    pose1_rdot1 = [pose1 deltaf[1]]
    pose2_rdot2 = [pose2 deltaf[2]]
    pose3_rdot3 = [pose3 deltaf[3]]
    pose4_rdot4 = [pose4 deltaf[4]]

    all_sats_scaled = vcat(pose1_rdot1,pose2_rdot2,pose3_rdot3,pose4_rdot4)

    if frequency[2] != 0

        pose5_rdot5 = [pose1 deltaf[5]]
        pose6_rdot6 = [pose2 deltaf[6]]
        pose7_rdot7 = [pose3 deltaf[7]]
        pose8_rdot8 = [pose4 deltaf[8]]

        all_sats_scaled = vcat(pose1_rdot1,pose2_rdot2,pose3_rdot3,pose4_rdot4, pose5_rdot5,pose6_rdot6,pose7_rdot7,pose8_rdot8)


    end

    #not even used 
    #Generate a guess at the centroid of all 4 sats on the surface of the Earth. Scaled to km
    centroid_guess = [(combined_eci[1,i]+combined_eci[7,i]+combined_eci[13,i]+combined_eci[19,i])/4, (combined_eci[2,i]+combined_eci[8,i]+combined_eci[14,i]+combined_eci[20,i])/4, (combined_eci[3,i]+combined_eci[9,i]+combined_eci[15,i]+combined_eci[21,i])/4] 
    onearth = sECEFtoGEOC(centroid_guess, use_degrees = true)
    geodetic = [onearth[1], onearth[2], 0]

    xyz_guess = sGEOCtoECEF(geodetic, use_degrees = true)*1e-3 #change to km
    
    mean_rescaled, all_r0, iters = tag_solve_Doppler(all_sats_scaled, xyz_guess, zenith_angles, time, zdot, frequency, d, r0_scaled)

    return mean_rescaled, all_r0, iters
end





