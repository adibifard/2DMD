# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 22:28:10 2022

@author: anshu
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy
           
sigma=0.34*10**(-9)
T=120
k=1.38*10**(-23)
E_k=k*T
N=108
m=6.7*10**(-26)
rho_s=0.8
rho=rho_s/(sigma**3)
V=N/rho_s
L=V**(1./3.)
l=L/(3)
M=3
print(L,l)

def position(L, N, M):
    ## Distribute N particles using fcc lattice over M*L cubic box 
    pos = np.zeros((N, 3),dtype=float)  # Init position array
    a = float(L)/M                      # Lattice constant
    j = 0                               # Set counter over all particles
    for ii in range(0, M):             # Unit cells, x-direction
        for iii in range(0, M):        # Unit cells, y-direction
            for iiii in range(0, M):   # Unit cells, z-direction
                ## Corner particle in unit cell
                pos[j,0] = a*ii         # x-position
                pos[j,1] = a*iii        # y-position
                pos[j,2] = a*iiii       # z-position
                j = j + 1               # Move to next particle
                ## First center particle in unit cell
                pos[j,0] = a*ii         
                pos[j,1] = a*iii+0.5*a   
                pos[j,2] = a*iiii+0.5*a  
                j = j + 1                
                ## Second center particle in unit cell
                pos[j,0] = a*ii+0.5*a    
                pos[j,1] = a*iii         
                pos[j,2] = a*iiii+0.5*a  
                j = j + 1                
                ## Third center particle in unit cell
                pos[j,0] = a*ii+0.5*a    
                pos[j,1] = a*iii+0.5*a  
                pos[j,2] = a*iiii       
                j = j + 1
    #pos = np.add(a*0.05,pos)
    print(j)
    return pos;

print(position(L,N,M))

## Initial velocity function
def vel(N):
	init_velocity = np.zeros((N,3),dtype=float)
	## Assign velocity from uniform random number
	for ii in range(0,N):
		init_velocity[ii,0] = np.random.normal(0.0,1.0)
		init_velocity[ii,1] = np.random.normal(0.0,1.0)
		init_velocity[ii,2] = np.random.normal(0.0,1.0)
	return init_velocity;

print(vel(N))

def normalize_vel(N,vel,T,E_k=-1):
    total_vel = np.zeros((N,3),dtype=float)
    ## Normalize velocities
    total_vel[:,0] = sum(vel[:,0])                    # total momentum in the x direction
    total_vel[:,1] = sum(vel[:,1])                    # total momentum in the y direction
    total_vel[:,2] = sum(vel[:,2])                    # total momentum in the z direction
    vel[:,0] = vel[:,0]-total_vel[:,0]/N         # readjust momentum per particle in the x direction
    vel[:,1] = vel[:,1]-total_vel[:,1]/N         # readjust momentum per particle in the y direction
    vel[:,2] = vel[:,2]-total_vel[:,2]/N         # readjust momentum per particle in the z direction
    print("sum=",sum(vel[:,0]),sum(vel[:,1]),sum(vel[:,2]))
    ## Determine rescaling factor
    if E_k == -1:       # if no total kinetic energy is given
        rescaling_factor = np.sqrt((3*(N-1)*T)/(sum(sum(np.array(vel)**2))))
    else:               # if total kinetic energy is already determined
        rescaling_factor = np.sqrt((3*(N-1)*T)/(2*E_k))
    vel = rescaling_factor*vel

    return vel;

print(normalize_vel(N,vel(N),T,E_k))

def pair_distance(Pos1,Pos2,L):
	distance = Pos2-Pos1 	#Distance in the x-dimension between position 1 and 2. Does not yet take into account periodic boundary conditions
	if np.fabs(distance)>(L/2):		# Check if there is a shorter distance with periodic boundary conditions
		if distance>0:				# Check for the sign of the distance
			distance = distance - L # and compensate accordingly
		else:
			distance = distance + L
	return distance;

def acceleration(N,pos,L):
	## Reset variables
    acceleration = np.zeros((N, 3),dtype=float)
    potential = np.zeros((N, 1),dtype=float)
    state = np.zeros((N,1),dtype=float)
    distance = np.zeros((3),dtype=float)
    dist_list = np.zeros((N,N),dtype=float)
	#virial = 0
    for ii in range(0, N):								# Loop over each particle
        for iii in range(0, N):						# Loop over all other particles
            if ii != iii:								# Exclude same particles
                distance[0] = pair_distance(pos[ii,0],pos[iii,0],L) 	# x distance
                distance[1] = pair_distance(pos[ii,1],pos[iii,1],L) 	# y distance
                distance[2] = pair_distance(pos[ii,2],pos[iii,2],L) 	# z distance
                abs_distance = np.sqrt(np.power(distance[0],2)+np.power(distance[1],2)+np.power(distance[2],2)) # Total distance
                dist_list[ii,iii] = abs_distance																# Save distance for correlation length calc
                V = 4*(np.power(abs_distance,-12)-np.power(abs_distance,-6))# Calculate lennard jones potential
                F = (48*np.power(abs_distance,-13)-24*np.power(abs_distance,-7))
                psi = (24*np.power(abs_distance,-6)*(2*np.power(abs_distance,-6)-1))
# 				if abs_distance < cutoff:																		# Set forces to zero if the distance is greater than the cutoff
# 					V = 4*(np.power(abs_distance,-12)-np.power(abs_distance,-6))			# Calculate lennard jones potential
# 					F = (48*np.power(abs_distance,-13)-24*np.power(abs_distance,-7))
# 				else:
# 					F = 0
# 					V = 0
# 				virial= abs_distance*F*0.5+virial										# Virial is the distance between two particles times the force between them, divide by 2 to avoid double counting
                acceleration[ii,0] = -F*(distance[0]/abs_distance)+acceleration[ii,0]	# Multiply the abs force with each component of the distance between two particles
                acceleration[ii,1] = -F*(distance[1]/abs_distance)+acceleration[ii,1]
                acceleration[ii,2] = -F*(distance[2]/abs_distance)+acceleration[ii,2]
            if ii < iii :
                potential[ii] = V + potential[ii]
                state[ii] = psi + state[ii]

	#dist_list = dist_list.flatten()							# Change the 2d array to a 1d array
	#dist_list = dist_list[dist_list !=0]					# remove zeroes from the 1d array
	#dist_hist,_ = np.histogram(dist_list,bins=hist_bins)	# histogram of distances, with width delta r.
    return acceleration*E_k/(m),potential*E_k;

def velocity_verlet(N, h, pos, v_0, a_0, L):
    v_half_h = np.add(v_0,0.5*a_0*h)
    pos_h = np.add(pos,v_half_h*h)
    
    # Impose periodic boundary condition
    pos_h = np.where(pos_h<0,L+pos_h,pos_h)             # If particle has a position smaller than 0, add L
    pos_h = np.where(pos_h>L, pos_h-L, pos_h)           # If particle has a position larger than L, subtract L

    a_h,potential = acceleration(N, pos_h, L)
    v_h = np.add(v_half_h,(0.5*a_h*h))

    # Calculate diffusion constant from the displacement
    #D = diffusion_constant(v_half_h,h,L)

    # Update acceleration, position and velocity to new values
    a_0 = a_h
    pos = pos_h
    v_0 = v_h

    return pos,v_0,a_0,potential;

pos = position( L,N,M )            # position
velocity = vel(N)       # momentum
a_0 = np.zeros((N,3),dtype=float)       # acceleration
time_dur=2000
h=0.01
time_step=np.zeros((time_dur),dtype=float)
kin_energy=np.zeros((time_dur),dtype=float)
pot_energy=np.zeros((time_dur),dtype=float)
tot_energy=np.zeros((time_dur),dtype=float)
## Time evolution
for t in range(0, time_dur):
    
    ## Equilibration phase
    time_step[t] = t
    ## Velocity verlet
    pos,velocity,a_0,potential = velocity_verlet( N, h, pos, velocity, a_0, L)
    ## Kinetic energy
    pot_energy[t] = sum(potential)
    kin_energy[t] = sum(sum(0.5*(np.power(velocity,2))))
    tot_energy[t] = kin_energy[t]+pot_energy[t]
    ## Intantanious temperature
    # T[t] = (float(2)/(3*(N-1)))*float(kin_energy[t])
    # ## Normalize momentum (with rescaling)
    # if np.mod(t,40) == 0 and t<=1201:
    #     velocity = normalize_momentum(N,velocity,T_d,kin_energy[t])
        
    # ## Equilibrium phase
    # if t>t_equil:
    #     pot_energy[t] = (0.5*sum(potential))/N                  # Potential energy per particle. Factor 0.5 to avoid double counting
    #     total_energy[t] = np.add(kin_energy[t],N*pot_energy[t])   # Total energy in the system
    #     P[t] = (phys.pressure(T[t],N,L,virial,r_c))/(T[t]*rho)  # Pressure
    #     sp_heat[t] = phys.specific_heat(N,T[t],kin_energy[t-50:t])     # Specific heat
    #     correlation_function = np.divide( ((2*np.power(L,3))/(N*(N-1)))*(np.mean(dist_hist,axis=1))/(4*np.pi*delta_r),np.power(np.multiply(hist_bins[1:],0.5),2))   # correlation function
        
    # ## Simulation progress
    # if np.mod(t,time_dur/100) == 0:
    #     # print ('%d%%' % t_prog)
    #     sys.stdout.write("Progress: %d%%   \r" % (t_prog) )
    #     sys.stdout.flush()
    #     t_prog = t_prog+1



## Plot data
#plot_data == 'y':
plt.plot(time_step,tot_energy)                          # Scale factor of 5.32E-8 to revert back to standard non-reduced units, m^2/s
plt.show()
