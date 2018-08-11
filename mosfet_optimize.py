#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 21:51:55 2018

@author: markditsworth
"""

import numpy as np
import subprocess
import ATLAS

def IV_data(filename):
    with open(filename,'rb') as fObj:
        data = fObj.readlines()
    
    data = data[20:]
    N = len(data)
    Id = np.zeros(N)
    Vg = np.zeros(N)
    Vd = np.zeros(N)
    
    for idx in range(N):
        line = data[idx].strip().split(' ')
        Id[idx] = float(line[9])
        Vg[idx] = float(line[4])
        Vd[idx] = float(line[7])
    
    return Id, Vg, Vd

# Sum of squares cost function
def cost(IdVdfile,IdVgfile):
    VTH_LIN_REF = 0
    VTH_TRANS_REF = 0
    VBD_REF = 0
    REF = np.array([VTH_LIN_REF,VTH_TRANS_REF,VBD_REF])
    # input and interpret log file
    vth_linear = 0
    vth_trans = 0
    vbd = 0
    output = np.array([vth_linear,vth_trans,vbd])
    
    cost = np.sum(np.power(output-REF,2))
    return cost

# Generates semiconductor geometry
def buildDeck(particle,IdVg,IdVd,particle_num):
    n_width = particle[0]
    p_width = particle[1]
    n_drift_width = particle[2]
    n_plus_doping = particle[3]
    n_drift_doping = particle[4]
    dit = particle[5]
    
    TOTAL_WIDTH = 600
    TOTAL_HEIGHT = 6
    DEPTH = 2000
    
    NSUB_HEIGHT = 3
    
    P_HEIGHT = 0.3
    
    N_HEIGHT = 0.15
    
    Cox = 422.018e-12 / 2
    p_doping = 2e17
    es = 9.66 *8.85e-12
    p_source_width = TOTAL_WIDTH - n_width - p_width - n_drift_width
    source_contact_width= p_source_width + (n_width/2.)
    tox = (es * (p_source_width * DEPTH / 1e12) / Cox) * 1e6 #* 0.6
    
    boundary = TOTAL_HEIGHT - NSUB_HEIGHT
    
    nsub = "0,%f %f,%f %f,%f 0,%f"%(TOTAL_HEIGHT, TOTAL_WIDTH,TOTAL_HEIGHT, TOTAL_WIDTH,TOTAL_HEIGHT-NSUB_HEIGHT, TOTAL_HEIGHT-NSUB_HEIGHT)
    
    ndrift = "%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f"%(0,TOTAL_HEIGHT-NSUB_HEIGHT, TOTAL_WIDTH,TOTAL_HEIGHT-NSUB_HEIGHT, TOTAL_WIDTH,0, TOTAL_WIDTH-n_drift_width,0, TOTAL_WIDTH-n_drift_width,P_HEIGHT, 0,P_HEIGHT)
    
    p = "%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f"%(0,P_HEIGHT, TOTAL_WIDTH-n_drift_width,P_HEIGHT, TOTAL_WIDTH-n_drift_width,0, TOTAL_WIDTH-n_drift_width-p_width,0, TOTAL_WIDTH-n_drift_width-p_width,N_HEIGHT, p_source_width,N_HEIGHT, p_source_width,0, 0,0)
    
    nsource = "%f,%f %f,%f %f,%f %f,%f %f,%f"%(p_source_width,N_HEIGHT, p_source_width+n_width,N_HEIGHT, p_source_width+n_width,0, source_contact_width,0, p_source_width,0)
    
    source = "%f,%f %f,%f %f,%f %f,%f %f,%f"%(0,0, source_contact_width,0, source_contact_width,-tox, source_contact_width,-0.3, 0,-0.3)
    
    oxide = "%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f"%(source_contact_width,0, p_source_width+n_width,0, TOTAL_WIDTH-n_drift_width,0, TOTAL_WIDTH,0, TOTAL_WIDTH,-tox, p_source_width+n_width-1,-tox, source_contact_width,-tox)
    
    gate = "%f,%f %f,%f %f,%f %f,%f"%(p_source_width+n_width-1,-tox, TOTAL_WIDTH,-tox, TOTAL_WIDTH,-tox-0.3, p_source_width+n_width-1,-tox-0.3)
    
    filename = "SiC_particle_%d.in"%particle_num
    
    ATLAS.deck(nsub,ndrift,p,nsource,source,oxide,gate,p_doping,n_drift_doping,n_plus_doping,dit,IdVg,IdVd,filename,boundary,particle_num,tox)
    
    return filename
    
def testIt():
    part=0
    for n_d in [1e15,5e15,1e16,5e16]:
        for n_p in [1e17,5e17,1e18,5e18]:
            fname = '%f_%f'%(n_d,n_p)
            p = [150,75,300,n_p,n_d,3e10]
            n = buildDeck(p,'IdVg_%s.log'%fname,'IdVd_%s.log'%fname,part)
            part += 1
            print n
        

    
def velocityUpdate(velocities,particle_vectors, particle_bests_v, global_best_v,r,s):
    # make global_best_V a matrix instead of a vector
    global_best_v = np.repeat(global_best_v.reshape((len(global_best_v),1)),particle_vectors.shape[1],axis=1)
    
    INERTIAL_CONST = 1 #may need to make a vector to account for scaling
    SOCIAL_COMP = 0.8
    COGNITIVE_COMP = 0.4
    
    # New velocity with inertial scaling
    vel = np.multiply(velocities,INERTIAL_CONST)
    # Add congnitive component to the velocity
    vel = vel + np.multiply(particle_bests_v-particle_vectors,r*COGNITIVE_COMP)
    # Add social component to the velocity
    vel = vel + np.multiply(global_best_v-particle_vectors,s*SOCIAL_COMP)
    return vel


def simulate(particles,run_num):
    particle_costs = np.zeros(particles.shape[1])
    
    for particle_num in range(particles.shape[1]):
        particle = particles[:,particle_num].flatten()
        IdVdfile = 'SiC_IdVd_run_%d_particle_%d.log'%(run_num,particle_num)
        IdVgfile = 'SiC_IdVg_run_%d_particle_%d.log'%(run_num,particle_num)
        # build deck commands
        deckFile = buildDeck(particle,IdVgfile,IdVdfile,particle_num)
        
        # call simulation
        cmd = '\\sedatools\\exe\\deckbuild -run %s'%deckFile
        subprocess.call(cmd.split(' '))
        
        # calculate cost
        particle_costs[particle_num] = cost(IdVdfile,IdVgfile)
    
    return particle_costs

def PSO():
    # CONSTANTS
    NUMBER_OF_PARTICLES = 10
    
    NUMBER_OF_ELEMENTS = 6
    
    # particle structure
    #
    # | n+width | Lch | N-drift/gate length | n+ | n- | dit | (in column vector form)
    ##########################################################
    
    particle = np.zeros((NUMBER_OF_ELEMENTS,NUMBER_OF_PARTICLES),dtype=np.float64)
    
    # initialize particles
    particle[0,0] = 180
    particle[1,0] = 190
    particle[2,0] = 200
    particle[3,0] = 1e17
    particle[4,0] = 1e16
    particle[5,0] = 3e10
    
    for x in range(NUMBER_OF_PARTICLES):
        if x>0:
            noise = 0.5 + np.random.random(size=(6,1))
            particle[:,x] = np.multiply(particle[:,0],noise)
    
    # initialize velocities
    v = np.random.rand(NUMBER_OF_ELEMENTS,NUMBER_OF_PARTICLES)*2 - 1
    
    # run initial simulations
    particle_costs = simulate(particle,0)
    
    # save particle-bests (intially equivalent to the first round)
    particle_best_costs = particle_costs.copy()
    particle_best_vectors = particle.copy()
    
    # save global-best
    global_best_cost = np.min(particle_best_costs)
    global_best_vector = particle[:,np.argmin(particle_best_costs)]
    
    costs = [global_best_cost]
    # main loop
    for q in range(10):
        print 'Iteration %d:'%(q+1)
        
        # randomize r and s vectors
        r = np.random.rand(NUMBER_OF_ELEMENTS,NUMBER_OF_PARTICLES)
        s = np.random.rand(NUMBER_OF_ELEMENTS,NUMBER_OF_PARTICLES)
        
        # update velocities
        v = velocityUpdate(v, particle, particle_best_vectors, global_best_vector,r,s)
        
        # update particles
        particle = particle + v
        
        # run simulation for each particle
        particle_costs = simulate(particle,q+1)
        
        # save particle bests
        for i,particle_cost in enumerate(particle_costs):
            if particle_cost < particle_best_costs[i]:
                particle_best_costs[i] = particle_cost
                particle_best_vectors[:,i] = particle[:,i]
        
        # save global best
        global_best_cost = np.min(particle_best_costs)
        costs.append(global_best_cost)
        global_best_vector = particle_best_vectors[:,np.argmin(particle_best_costs)]
    