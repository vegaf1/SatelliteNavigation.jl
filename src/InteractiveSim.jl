#Functions for Interactive Simulation

import SatelliteNavigation as SN
#using GLMakie, FileIO, Colors
using Distributions


#Initialize Constants
const G = 6.67430e-20 #km**3/(kg * s**2)
const m_1 = 5.97219e24 #mass of the Earth in kg
const m_2 = 1 #mass of satellite in kg
const μ = G*m_1 #gravitational parameter
const R_earth = 6378.12 # radius of Earth in km
const h = 1 #time step in seconds
const hi = 350e3 #in meters. Altitude of ionosphere thin shell approx.

mu = 8e-6
sigma = 3e-6
lb = 3e-6
ub = 13e-6
d = Truncated(Normal(mu,sigma), lb, ub)

#Define structures for TOA and Doppler Units

struct TOA_units

    distance_scale::Float64
    time_scale::Float64
    c::Float64

    function TOA_units()
        #TOA scales
        distance_scale = R_EARTH*(1e-3) # scales km to custom scale
        time_scale = 1/(C_LIGHT/R_EARTH) #scales seconds to custom scale
        c = 1.0 # speed of light
        new(distance_scale, time_scale, c)
    end

end

struct Doppler_units

    distance_scale::Float64
    time_scale::Float64
    velocity_scale::Float64
    c::Float64

    function Doppler_units()

        #Doppler scales
        distance_scale = R_EARTH*(1e-3) #scales km to custom scale
        time_scale = (1/(C_LIGHT/R_EARTH))*1e4 #scales seconds to custom scale
        velocity_scale = distance_scale/time_scale #scales km/s to custom units
        c = 1.0*1e4 #speed of light
        new(distance_scale, time_scale, velocity_scale, c)

    end

end

function sat_dynamics(x)
    
    sat_pose = SA[x[1], x[2], x[3]]
    
    #This WORKS
    upperleft = @SMatrix zeros(3,3)
    upperright = SMatrix{3,3}(1.0I)
    lowerleft = SMatrix{3,3}(1I)*(-μ/(norm(sat_pose))^3)
    lowerright = @SMatrix zeros(3,3)
    
    
    upper = [upperleft upperright]
    lower = [lowerleft lowerright]
    
    A = [upper; lower]
    
    xdot = A*x
    
    return xdot
    
end

#Runga Kutta 4th order algorithm to find the next state

function RK4_orbit(x, h)
    
    f1_ = sat_dynamics(x)
    f2_ = sat_dynamics(x+0.5*h*f1_)
    f3_ = sat_dynamics(x+0.5*h*f2_)
    f4_ = sat_dynamics(x+h*f3_)
    
    xnext = x+(h/6.0)*(f1_+2*f2_+2*f3_+f4_)

    return xnext
    
end

#Propagate the orbit

function twobody_rk4(fun, x0, Tf, h)
    
    t = Array(range(0,Tf, step=h)) #create a range to final time at constant time step
    
    all_x = zeros(length(x0), length(t)) #variable to store all x
    
    all_x[:,1] .= x0 #set the initial state
    
    for k=1:(length(t) - 1)
        all_x[:,k+1] .= RK4_orbit(all_x[:,k], h)
    end
    
    return all_x, t
end

function tag_xyz(latitude, toa_units, doppler_units)

    #Setting a known tag position
    #Fix the longitude and altitude, vary the latitude

    #longitude, latitude, altitude
    tag_geof = [-165.4545, latitude, 0]

    #gives tag position in meters
    tag = sGEOCtoECEF(tag_geof, use_degrees = true)

    #Equitorial Position
    #Rescale to km
    x0= tag[1]*1e-3
    y0 = tag[2]*1e-3
    z0 = tag[3]*1e-3
    t0 = 1e-5
    bdot = 1e4*1e-3 #convert to km/s almost working

    #TOA unknowns unscaled 
    r0 = [x0, y0, z0, t0]

    #TOA unknowns scaled
    r0_scaled = [x0/toa_units.distance_scale, y0/toa_units.distance_scale, z0/toa_units.distance_scale, t0/toa_units.time_scale]

    #knowns for Doppler case
    TEC1 = 3e-4
    TEC2 = 4e-4
    TEC3 = 3.5e-4
    TEC4 = 4.5e-4

    #Doppler unknowns scaled
    r0_scaled_Doppler = [x0/doppler_units.distance_scale, y0/doppler_units.distance_scale, z0/doppler_units.distance_scale, bdot/doppler_units.velocity_scale, TEC1, TEC2, TEC3, TEC4]
    
    return r0, r0_scaled, r0_scaled_Doppler

end

#update the orbit when the slider changes

function orbit_update(trueanom, RAAN_sep, delta_sep, altitude)

    h=1 #timestep
    
    iss1 = [6371e3 + (altitude*1e3), 0.0004879, 90.6391, 194.0859- (RAAN_sep/2), 151.2014, 190];
    iss2 = [6371e3 + (altitude*1e3), 0.0004879, 90.6391, 194.0859 - (RAAN_sep/2), 151.2014, 190+trueanom];
    iss3 = [6371e3 + (altitude*1e3), 0.0004879, 90.6391, 195.0859 + (RAAN_sep/2), 151.2014, 190+delta_sep]; 
    iss4 = [6371e3 + (altitude*1e3), 0.0004879, 90.6391, 195.0859 + (RAAN_sep/2), 151.2014, 190+delta_sep+trueanom]; 
    
    eci0_1 = sOSCtoCART(iss1, use_degrees=true)
    eci0_2 = sOSCtoCART(iss2, use_degrees=true)
    eci0_3 = sOSCtoCART(iss3, use_degrees=true)
    eci0_4 = sOSCtoCART(iss4, use_degrees=true)
    
    T = orbit_period(iss1[1])
    
    #in km
    x0_1 = SA[eci0_1[1], eci0_1[2], eci0_1[3],eci0_1[4], eci0_1[5], eci0_1[6]]*1e-3 
    x0_2 = SA[eci0_2[1], eci0_2[2], eci0_2[3],eci0_2[4], eci0_2[5], eci0_2[6]]*1e-3 
    x0_3 = SA[eci0_3[1], eci0_3[2], eci0_3[3],eci0_3[4], eci0_3[5], eci0_3[6]]*1e-3 
    x0_4 = SA[eci0_4[1], eci0_4[2], eci0_4[3],eci0_4[4], eci0_4[5], eci0_4[6]]*1e-3 
    
    #Find state of orbit 1,2,3,4
    eci_1, t_hist1 = twobody_rk4(sat_dynamics, x0_1, T, h)
    eci_2, t_hist2 = twobody_rk4(sat_dynamics, x0_2, T, h)
    eci_3, t_hist3 = twobody_rk4(sat_dynamics, x0_3, T, h)
    eci_4, t_hist4 = twobody_rk4(sat_dynamics, x0_4, T, h)
    
    #6557 is the size of the array when it is at 1200. Need to make it a fixed size to be able to make it 
    #interactive
    
    size_eci1 = size(eci_1)
    size_eci2 = size(eci_2)
    size_eci3 = size(eci_3)
    size_eci4 = size(eci_4)
    
    eci_1_sized = zeros(6,6557)
    eci_2_sized = zeros(6,6557)
    eci_3_sized = zeros(6,6557)
    eci_4_sized = zeros(6,6557)
    
    eci_1_sized[1:size_eci1[1], 1:size_eci1[2]] = eci_1
    eci_2_sized[1:size_eci2[1], 1:size_eci2[2]] = eci_2
    eci_3_sized[1:size_eci3[1], 1:size_eci3[2]] = eci_3
    eci_4_sized[1:size_eci4[1], 1:size_eci4[2]] = eci_4
    
    centroid_guess = [(eci_1[1,1]+eci_2[1,1]+eci_3[1,1]+eci_4[1,1])/4, (eci_1[2,1]+eci_2[2,1]+eci_3[2,1]+eci_4[2,1])/4, (eci_1[3,1]+eci_2[3,1]+eci_3[3,1]+eci_4[3,1])/4] 
    onearth = sECEFtoGEOC(centroid_guess, use_degrees = true)
    geodetic = [onearth[1], onearth[2], 0]

    #Guess
    xyz = sGEOCtoECEF(geodetic, use_degrees = true)*1e-3
    
    return eci_1_sized, eci_2_sized, eci_3_sized, eci_4_sized, xyz
        
end


function get_P(time_accuracy, frequency_accuracy, toa_units, doppler_units)
    #Errors in satellite pose
    #distance scale is the same for toa and doppler
    position_accuracy = 0.1 * 1e-3/toa_units.distance_scale #(0.1 meters) 2 sigma standard deviation
    velocity_accuracy = 0.01 *1e-3/doppler_units.velocity_scale #(0.01 meters/second) 2 sigma standard deviation
    
    #covariance matrix is made up by the squares of the standard deviation (variance)
    P_TOA_1sat = diagm([(position_accuracy)^2, (position_accuracy)^2, (position_accuracy)^2, (time_accuracy/toa_units.time_scale)^2])
    P_Doppler_1sat = diagm([(position_accuracy)^2, (position_accuracy)^2, (position_accuracy)^2, (velocity_accuracy)^2, (velocity_accuracy)^2, (velocity_accuracy)^2, (frequency_accuracy)^2])
    
    #Make the block diagonal for the 8 measurments
    
    P_TOA = BlockDiagonal([P_TOA_1sat, P_TOA_1sat, P_TOA_1sat, P_TOA_1sat, P_TOA_1sat, P_TOA_1sat, P_TOA_1sat, P_TOA_1sat])
    P_Doppler = BlockDiagonal([P_Doppler_1sat , P_Doppler_1sat , P_Doppler_1sat , P_Doppler_1sat , P_Doppler_1sat , P_Doppler_1sat , P_Doppler_1sat , P_Doppler_1sat ])

    return P_TOA, P_Doppler
end



function LS_Amatrix(x, y, scale, A, angle, toa_units)
    
    jacobian_pose = (x[1:3] - A*y[1:3])/norm(x[1:3]-A*y[1:3])
    
    dpdtau = toa_units.c
    dpdTEC = (scale/cosd(angle))*(1e-3)/toa_units.distance_scale
    
    return [jacobian_pose; dpdtau; dpdTEC]
end

function LS_Amatrix_dpdy(x, y, scale, A, angle, toa_units)
    
    #c = C_LIGHT
    
    #works
    jacobian_pose = (y[1:3] - A'*x[1:3])/norm(x[1:3]-A*y[1:3])
    
    dpdtau = -toa_units.c
    
    return [jacobian_pose; dpdtau]
end

function customjacobian(initial_x0,sat_pose_t, zenith_angles, frequency, toa_units)

    # c=1
    
    while norm(SN.TOA_residual(initial_x0, sat_pose_t, zenith_angles, frequency))>1e-10
        
        omega = OMEGA_EARTH*toa_units.time_scale

        angles = zeros(4)

        for i=1:4
            angles[i] = asind((R_EARTH*sind(zenith_angles[i]))/(R_EARTH + hi))
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
        
        P = SN.TOA_residual(initial_x0, sat_pose_t, zenith_angles, frequency)

        deltaP = -P

        #Create the A matrix to solve the LS problem
        
        A = zeros(8,8)

        A_row1 = LS_Amatrix(initial_x0, sat_pose_t[1,:], scale_f1, A1, angles[1], toa_units)'
        A_row2 = LS_Amatrix(initial_x0, sat_pose_t[2,:], scale_f1, A2, angles[2], toa_units)'
        A_row3 = LS_Amatrix(initial_x0, sat_pose_t[3,:], scale_f1, A3, angles[3], toa_units)'
        A_row4 = LS_Amatrix(initial_x0, sat_pose_t[4,:], scale_f1, A4, angles[4], toa_units)'
        A_row5 = LS_Amatrix(initial_x0, sat_pose_t[5,:], scale_f2, A5, angles[1], toa_units)'
        A_row6 = LS_Amatrix(initial_x0, sat_pose_t[6,:], scale_f2, A6, angles[2], toa_units)'
        A_row7 = LS_Amatrix(initial_x0, sat_pose_t[7,:], scale_f2, A7, angles[3], toa_units)'
        A_row8 = LS_Amatrix(initial_x0, sat_pose_t[8,:], scale_f2, A8, angles[4], toa_units)'

        
        A[1,:] = [A_row1[1:5]; zeros(3)]
        A[2,:] = [A_row2[1:4];0.0; A_row2[5]; zeros(2)]
        A[3,:] = [A_row3[1:4];zeros(2); A_row3[5]; 0.0]
        A[4,:] = [A_row4[1:4];zeros(3); A_row4[5]]
        A[5,:] = [A_row5[1:5]; zeros(3)]
        A[6,:] = [A_row6[1:4];0.0; A_row6[5]; zeros(2)]
        A[7,:] = [A_row7[1:4];zeros(2); A_row7[5]; 0.0]
        A[8,:] = [A_row8[1:4];zeros(3); A_row8[5]]
        
        F = qr(A)

        b = transpose(F.Q)*deltaP
        
        deltax = F.R\b
        
        initial_x0 = initial_x0+deltax

    end

    x_sol = initial_x0
    
    return x_sol
    
end

#Used to create the A matrix to propogate the covariance of the satellites to the covariance of the tag

function createA_TOA(sat_pose_t, zenith_angles, x_sol, frequency, toa_units)
    
    omega = OMEGA_EARTH*toa_units.time_scale #scaled

    angles = zeros(4)

    for i=1:4
        angles[i] = asind((R_EARTH*sind(zenith_angles[i]))/(R_EARTH + hi))
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

    #create dfdy
    dfdx = zeros(8,8)

    dfdx_row1 = LS_Amatrix(x_sol, sat_pose_t[1,:], scale_f1, A1, angles[1], toa_units)'
    dfdx_row2 = LS_Amatrix(x_sol, sat_pose_t[2,:], scale_f1, A2, angles[2], toa_units)'
    dfdx_row3 = LS_Amatrix(x_sol, sat_pose_t[3,:], scale_f1, A3, angles[3], toa_units)'
    dfdx_row4 = LS_Amatrix(x_sol, sat_pose_t[4,:], scale_f1, A4, angles[4], toa_units)'
    dfdx_row5 = LS_Amatrix(x_sol, sat_pose_t[5,:], scale_f2, A5, angles[1], toa_units)'
    dfdx_row6 = LS_Amatrix(x_sol, sat_pose_t[6,:], scale_f2, A6, angles[2], toa_units)'
    dfdx_row7 = LS_Amatrix(x_sol, sat_pose_t[7,:], scale_f2, A7, angles[3], toa_units)'
    dfdx_row8 = LS_Amatrix(x_sol, sat_pose_t[8,:], scale_f2, A8, angles[4], toa_units)'

    dfdx[1,:] = [dfdx_row1[1:5]; zeros(3)]
    dfdx[2,:] = [dfdx_row2[1:4];0.0; dfdx_row2[5]; zeros(2)]
    dfdx[3,:] = [dfdx_row3[1:4];zeros(2); dfdx_row3[5]; 0.0]
    dfdx[4,:] = [dfdx_row4[1:4];zeros(3); dfdx_row4[5]]
    dfdx[5,:] = [dfdx_row5[1:5]; zeros(3)]
    dfdx[6,:] = [dfdx_row6[1:4];0.0; dfdx_row6[5]; zeros(2)]
    dfdx[7,:] = [dfdx_row7[1:4];zeros(2); dfdx_row7[5]; 0.0]
    dfdx[8,:] = [dfdx_row8[1:4];zeros(3); dfdx_row8[5]]

    #derivative wrt measurements (satellite position and times)
    dfdy = zeros(8,32)

    #create dfdx
    dfdy[1, 1:4] = LS_Amatrix_dpdy(x_sol, sat_pose_t[1,:], scale_f1, A1, angles[1], toa_units)'[1:4]
    dfdy[2, 5:8] = LS_Amatrix_dpdy(x_sol, sat_pose_t[2,:], scale_f1, A2, angles[2], toa_units)'[1:4]
    dfdy[3, 9:12] = LS_Amatrix_dpdy(x_sol, sat_pose_t[3,:], scale_f1, A3, angles[3], toa_units)'[1:4]
    dfdy[4, 13:16] = LS_Amatrix_dpdy(x_sol, sat_pose_t[4,:], scale_f1, A4, angles[4], toa_units)'[1:4]
    dfdy[5, 17:20] = LS_Amatrix_dpdy(x_sol, sat_pose_t[5,:], scale_f2, A5, angles[1], toa_units)'[1:4]
    dfdy[6, 21:24] = LS_Amatrix_dpdy(x_sol, sat_pose_t[6,:], scale_f2, A6, angles[2], toa_units)'[1:4]
    dfdy[7, 25:28] = LS_Amatrix_dpdy(x_sol, sat_pose_t[7,:], scale_f2, A7, angles[3], toa_units)'[1:4]
    dfdy[8, 29:32] = LS_Amatrix_dpdy(x_sol, sat_pose_t[8,:], scale_f2, A8, angles[4], toa_units)'[1:4]

    A = -inv(dfdx)*dfdy
    
    return A
end

#calculate the covariance of the tag from the TOA method 

function calculate_covariance_TOA(sat_poses, zenith_angles, guess, P_TOA, r0_scaled, frequency, toa_units, latitude)
    TEC = rand(d,4) 
    r0_TEC= vec([r0_scaled TEC])
    sat_poses_scaled = sat_poses/toa_units.distance_scale #scale the data
    t = SN.get_time(r0_TEC, sat_poses_scaled, zenith_angles, frequency)
    sat_poses_2 = vcat(sat_poses_scaled, sat_poses_scaled)
    sat_pose_t = [sat_poses_2 t]
    
    #initial guess
    
    initial_x0 = [guess[1]/toa_units.distance_scale, guess[2]/toa_units.distance_scale, guess[3]/toa_units.distance_scale, 0.0001/toa_units.time_scale, 3e-6, 3.5e-6, 2e-6, 4e-6]

    #get the converged solution
    x_sol = customjacobian(initial_x0,sat_pose_t, zenith_angles, frequency, toa_units)
    
    if all(isfinite, x_sol) == true
    
        A = createA_TOA(sat_pose_t, zenith_angles, x_sol, frequency, toa_units)

        #convert the covariance of the satellites to the covariance of the tag

        covariance_tag = A*P_TOA*transpose(A)

        #Only get the pose covariance
        pose_covariance = covariance_tag[1:3, 1:3]

        #scale the pose covariance to meters
        pose_covariance_scaled = covariance_tag[1:3, 1:3]*(1e3*toa_units.distance_scale)^2

        #longitude and latitude of the tag position. 
        lambda = -165.4545 #longitude 
        phi =  latitude #latitude

        #Switch from xyz to East North Up coordinate system

        R = [-sind(lambda) cosd(lambda) 0; 
            -sind(phi)*cosd(lambda) -sind(phi)*sind(lambda) cosd(phi);
            cosd(phi)*cosd(lambda) cosd(phi)*sind(lambda) sind(phi)]

        #Rotate the covariance matrix
        covariance_matrix_scaled_rotated = R*pose_covariance_scaled*R'

        #Find the standard deviation of the tag position
        PDOP_LS = sqrt(covariance_tag[1,1] + covariance_tag[2,2] + covariance_tag[3,3])

        #This is the standard deviation of the tag postion scaled in meters and seconds
        PDOP_LS_scaled = PDOP_LS*toa_units.distance_scale*1e3
 
    else
        
        PDOP_LS = 0
        PDOP_LS_scaled = eyes(8)
        
    end
    
    return PDOP_LS_scaled, covariance_matrix_scaled_rotated
    
end

#sat_pose_f is a vector of all the satellite positions, velocities. Vector in order to use Forward Diff

function doppler_residual_vec(x, sat_pose_f, time, z, zdot, frequency, doppler_units)
    
    omega = OMEGA_EARTH*doppler_units.time_scale #transform to custom scale
    
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
    
    #bdot is a velocity. therefore it is scaled
    #Ionosphere Errors at frequency 1
    Idot1 = (((40.3*x[5])*((R_EARTH*sind(z[1]))*(R_EARTH*cosd(z[1])*zdot[1])/(R_EARTH + hi)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[1]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale
    Idot2 = (((40.3*x[6])*((R_EARTH*sind(z[2]))*(R_EARTH*cosd(z[2])*zdot[2])/(R_EARTH + hi)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[2]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale
    Idot3 = (((40.3*x[7])*((R_EARTH*sind(z[3]))*(R_EARTH*cosd(z[3])*zdot[3])/(R_EARTH + hi)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[3]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale
    Idot4 = (((40.3*x[8])*((R_EARTH*sind(z[4]))*(R_EARTH*cosd(z[4])*zdot[4])/(R_EARTH + hi)^2))/((frequency[1]^2)*(1-((R_EARTH*sind(z[4]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale

    #Ionosphere Errors at frequency 2
    Idot5 = (((40.3*x[5])*((R_EARTH*sind(z[1]))*(R_EARTH*cosd(z[1])*zdot[1])/(R_EARTH + hi)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[1]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale
    Idot6 = (((40.3*x[6])*((R_EARTH*sind(z[2]))*(R_EARTH*cosd(z[2])*zdot[2])/(R_EARTH + hi)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[2]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale
    Idot7 = (((40.3*x[7])*((R_EARTH*sind(z[3]))*(R_EARTH*cosd(z[3])*zdot[3])/(R_EARTH + hi)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[3]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale
    Idot8 = (((40.3*x[8])*((R_EARTH*sind(z[4]))*(R_EARTH*cosd(z[4])*zdot[4])/(R_EARTH + hi)^2))/((frequency[2]^2)*(1-((R_EARTH*sind(z[4]))/(R_EARTH+hi))^2)^1.5))*1e20*1e-3/doppler_units.velocity_scale


    #Measurments at frequency 1
    res1 = (frequency[1]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A1*sat_pose_f[1:3]))*(-2*x[1:3]'*A1*omegahat*sat_pose_f[1:3] - 2*x[1:3]'*A1*sat_pose_f[4:6]+sat_pose_f[1:3]'*sat_pose_f[4:6]+sat_pose_f[4:6]'*sat_pose_f[1:3]) + Idot1 + x[4]) - sat_pose_f[7]
    res2 = (frequency[1]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A2*sat_pose_f[8:10]))*(-2*x[1:3]'*A2*omegahat*sat_pose_f[8:10] - 2*x[1:3]'*A2*sat_pose_f[11:13]+sat_pose_f[8:10]'*sat_pose_f[11:13]+sat_pose_f[11:13]'*sat_pose_f[8:10]) + Idot2 + x[4]) - sat_pose_f[14]
    res3 = (frequency[1]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A3*sat_pose_f[15:17]))*(-2*x[1:3]'*A3*omegahat*sat_pose_f[15:17] - 2*x[1:3]'*A3*sat_pose_f[18:20]+sat_pose_f[15:17]'*sat_pose_f[18:20]+sat_pose_f[18:20]'*sat_pose_f[15:17]) + Idot3 + x[4]) - sat_pose_f[21]
    res4 = (frequency[1]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A4*sat_pose_f[22:24]))*(-2*x[1:3]'*A4*omegahat*sat_pose_f[22:24] - 2*x[1:3]'*A4*sat_pose_f[25:27]+sat_pose_f[22:24]'*sat_pose_f[25:27]+sat_pose_f[25:27]'*sat_pose_f[22:24]) + Idot4 + x[4]) - sat_pose_f[28]
    
    #Measurments at frequency 2
    res5 = (frequency[2]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A1*sat_pose_f[1:3]))*(-2*x[1:3]'*A1*omegahat*sat_pose_f[1:3] - 2*x[1:3]'*A1*sat_pose_f[4:6]+sat_pose_f[1:3]'*sat_pose_f[4:6]+sat_pose_f[4:6]'*sat_pose_f[1:3]) + Idot5 + x[4]) - sat_pose_f[35]
    res6 = (frequency[2]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A2*sat_pose_f[8:10]))*(-2*x[1:3]'*A2*omegahat*sat_pose_f[8:10] - 2*x[1:3]'*A2*sat_pose_f[11:13]+sat_pose_f[8:10]'*sat_pose_f[11:13]+sat_pose_f[11:13]'*sat_pose_f[8:10]) + Idot6 + x[4]) - sat_pose_f[42]
    res7 = (frequency[2]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A3*sat_pose_f[15:17]))*(-2*x[1:3]'*A3*omegahat*sat_pose_f[15:17] - 2*x[1:3]'*A3*sat_pose_f[18:20]+sat_pose_f[15:17]'*sat_pose_f[18:20]+sat_pose_f[18:20]'*sat_pose_f[15:17]) + Idot7 + x[4]) - sat_pose_f[49]
    res8 = (frequency[2]/doppler_units.c)*(0.5*(1/norm(x[1:3] - A4*sat_pose_f[22:24]))*(-2*x[1:3]'*A4*omegahat*sat_pose_f[22:24] - 2*x[1:3]'*A4*sat_pose_f[25:27]+sat_pose_f[22:24]'*sat_pose_f[25:27]+sat_pose_f[25:27]'*sat_pose_f[22:24]) + Idot8 + x[4]) - sat_pose_f[56]

    
    return [res1; res2; res3; res4; res5; res6; res7; res8]
   
end

#Create the A matrix for Doppler to obtain the covariance of the tag

function createA_Doppler(sat_pose_t, zenith_angles, x_sol, time, zdot, frequency, doppler_units)
    satpose_vec = vec(sat_pose_t')
    dfdy_fd = ForwardDiff.jacobian(dall_sats_scaled -> doppler_residual_vec(x_sol, dall_sats_scaled, time, zenith_angles, zdot, frequency, doppler_units), satpose_vec)
    dfdx_fd = ForwardDiff.jacobian(dx -> doppler_residual_vec(dx, satpose_vec, time, zenith_angles, zdot, frequency, doppler_units), x_sol)
    
    #use implicit function theorem
    A = -inv(dfdx_fd)*dfdy_fd
    
    return A
end

#Solve the Doppler LS problem
function customjacobian_Doppler(initial_x0, sat_pose_t, time, zenith_angles, zdot, frequency, doppler_units)

    #make sure sat_pose_t is a vector
    satpose_vec = vec(sat_pose_t')
    
    while norm(doppler_residual_vec(initial_x0, satpose_vec, time, zenith_angles, zdot, frequency, doppler_units))>1e-10
        
        res = doppler_residual_vec(initial_x0, satpose_vec, time, zenith_angles, zdot, frequency, doppler_units)
    
        deltaP = -res
        
        A = ForwardDiff.jacobian(dx -> doppler_residual_vec(dx, satpose_vec, time, zenith_angles, zdot, frequency, doppler_units), initial_x0)
        
        F = qr(A)

        d_ = transpose(F.Q)*deltaP

        deltax = F.R\d_
        
        initial_x0 = initial_x0+deltax
        
    end
    
    x_sol = initial_x0
    
    return x_sol
        
end

function calculate_covariance_Doppler(pose1, pose2, pose3, pose4, zenith_angles, zdot, P_Doppler, r0_scaled_Doppler, frequency, doppler_units, latitude)
    
    #the difference added to actual solution to converge (for Doppler scenerio)
    diff = 50e3*1e-3/doppler_units.distance_scale 
    
    #assume a propogation time to include rotation of earth
    time = [0.006, 0.006, 0.006, 0.008]/doppler_units.time_scale
    
    #initial guess for Doppler algorithm
    initial_x0 = [r0_scaled_Doppler[1]+diff, r0_scaled_Doppler[2]+diff, r0_scaled_Doppler[3]+diff, 1.002*1e4*1e-3/doppler_units.velocity_scale, 3.1e-4, 4.2e-4, 3.2e-4, 4.7e-4]
    
    sat_poses = [pose1; pose2; pose3; pose4]

    #Obtain the true measurments for the doppler case
    deltaf = SN.doppler_measurment(r0_scaled_Doppler, sat_poses, time, zenith_angles, zdot, frequency)
    
    #augment the frequency measurments to the satellite position and velocities
    pose1_rdot1 = [pose1 deltaf[1]]
    pose2_rdot2 = [pose2 deltaf[2]]
    pose3_rdot3 = [pose3 deltaf[3]]
    pose4_rdot4 = [pose4 deltaf[4]]
    
    pose5_rdot5= [pose1 deltaf[5]]
    pose6_rdot6= [pose2 deltaf[6]]
    pose7_rdot7= [pose3 deltaf[7]]
    pose8_rdot8= [pose4 deltaf[8]]
    
    all_sats_scaled = vcat(pose1_rdot1,pose2_rdot2,pose3_rdot3,pose4_rdot4, pose5_rdot5,pose6_rdot6,pose7_rdot7,pose8_rdot8)
    
    #Find the Least Squares solution
    x_sol = customjacobian_Doppler(initial_x0, all_sats_scaled, time, zenith_angles, zdot, frequency, doppler_units)
    
    if all(isfinite, x_sol) == true
        
    
        #Get the A matrix to find the covaraince of the tag
        A = createA_Doppler(all_sats_scaled, zenith_angles, x_sol, time, zdot, frequency, doppler_units)

        #Calculate the Covariance of the tag
        covariance_tag = A*P_Doppler*transpose(A)

        #Only get the part of the covariance matrix that has to do with position
        pose_covariance = covariance_tag[1:3, 1:3]

        #Scale back to meters and seconds
        pose_covariance_scaled = covariance_tag[1:3, 1:3]*(1e3*doppler_units.distance_scale)^2

        lambda = -165.4545 #longitude 
        phi =  latitude #latitude

        #Rotation matrix to update the covariance from xyz to east north up coordinate system
        R = [-sind(lambda) cosd(lambda) 0; 
            -sind(phi)*cosd(lambda) -sind(phi)*sind(lambda) cosd(phi);
            cosd(phi)*cosd(lambda) cosd(phi)*sind(lambda) sind(phi)]

        covariance_matrix_scaled_rotated = R*pose_covariance_scaled*R'

        PDOP_LS = sqrt(covariance_tag[1,1] + covariance_tag[2,2] + covariance_tag[3,3])

        PDOP_LS_scaled = PDOP_LS*doppler_units.distance_scale*1e3

        #return PDOP_LS_scaled, covariance_matrix_scaled_rotated
        
    else
        
        PDOP_LS_scaled = 0
        covariance_matrix_scaled_rotated = eyes(8)
    end
        
        
    return PDOP_LS_scaled, covariance_matrix_scaled_rotated
    
end

#DRAW THE ELLIPSE
function draw_ellipse(covariance_matrix)
    
    eig_value = eigvals(covariance_matrix)
    eig_vecs = eigvecs(covariance_matrix)
    
    max_eigval_index = argmax(eig_value)
    min_eigval_index = argmin(eig_value)

    max_eigval = eig_value[max_eigval_index]
    min_eigval = eig_value[min_eigval_index]
    
    max_eigvec = eig_vecs[:, max_eigval_index]
    min_eigvec = eig_vecs[:, min_eigval_index]
    
    angle = atan(max_eigvec[2], max_eigvec[1])
    
    if angle<0
        angle = angle + 2*pi;
    end
    
    #95% confidence interval
    chisquare_val = 2.4477
    theta_grid = 0:0.01:2*pi

    phi = angle;

    a = chisquare_val*sqrt(max_eigval)
    b = chisquare_val*sqrt(min_eigval)

    ellipse_x = zeros(size(theta_grid)[1], 1)
    ellipse_y = zeros(size(theta_grid)[1], 1)

    for i=1:size(theta_grid)[1]

        #Get ellipse coordinates 
        ellipse_x[i,1] = a*cos(theta_grid[i])
        ellipse_y[i,1] = b*sin(theta_grid[i])

    end
    
    R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
    r_ellipse = [ellipse_x ellipse_y]*R
    
    return r_ellipse
    
end



















