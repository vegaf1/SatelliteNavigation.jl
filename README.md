# SatelliteNavigation.jl

## Background
A Low-Earth Orbiting (LEO) satellite constellation is simulated to track a target on the surface of the Earth using two distinct methods. The on-Earth target is assumed to have a transmitting tag that communicates with the satellites at a certain frequency to obtain ranging measurements.   The first simulation method utilizes time of flight measurements to measure the relative range between each satellite in the constellation and the transmitting tag. The second method uses the Doppler shift which uses the change between the transmitting and receiving frequency to calculate the distance between the satellite and the tag. Errors in the satellite ephemeris data, clock, and frequency measurements are taken into account in the simulation by representing them as random variables which are sampled from a scaled Gaussian distribution. Satellites with multiple antennas were also investigated using a time difference of arrival method; however, it achieved poor position accuracy due to the high uncertainty in cubesat attitude. The result of the simulation provided the uncertainty in the tag position in the form of a covariance ellipse. A change in orbital parameters and constellation design heavily influenced the tag estimate uncertainty. The change in pseudorange due to the ionosphere was also considered in the simulation, as this effect can affect the overall estimate by up to 40m depending on the time of day. 

## Orbit Design
The following parameters were changed in the simulation: 
- RAAN seperation (\omega_{sep})
- True anomaly seperation (\theta_{sep})
- Delta True Anomaly seperation (\Delta \theta_{sep})
- Orbit Altitude 

![Alt text](/home/fausto/satellite_formation.png "Satellite Constellation")

