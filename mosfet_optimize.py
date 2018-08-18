#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 21:51:55 2018

@author: markditsworth
"""
# Make Model functions belong to a class to better integration of temperature
import numpy as np
import subprocess
import ATLAS
from scipy.optimize import curve_fit

def linear(x,m,b):
    y = m*x + b
    return y

def deriv(x,y,window=10):
    x_ = []
    dy = []
    
    assert len(x) == len(y), 'Data mismatch. Must be of same length!'
    
    for i in range(len(x)):
        if i >= (window/2) and i < (len(x)-(window/2)):
            x_.append(x[i])
            
            yy = y[i-(window/2):i+((window/2)-1)]
            xx = x[i-(window/2):i+((window/2)-1)]
            
            popt,_ = curve_fit(linear,xx,yy)
            
            dy.append(popt[0])
        
    return np.array(x_), np.array(dy)

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

def vth_linear(i_,v_,window_size=10):
    v = v_[v_>3]
    i = i_[v_>3]
    sample_idx = np.arange(window_size)
    lines = []
    Rs = []
    while sample_idx[-1] < len(v):
        sample_i = i[sample_idx]
        sample_v = v[sample_idx]
        popt,_ = curve_fit(linear,sample_v,sample_i)
        lines.append(popt)
        # calculate r-squared
        residuals = sample_i- linear(sample_v, popt[0],popt[1])
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((sample_i-np.mean(sample_i))**2)
        Rs.append( 1 - (ss_res / ss_tot))
        sample_idx = sample_idx + 1
    best_line = lines[np.argmax(np.array(Rs))]
    
    m = best_line[0]
    b = best_line[1]
    vth = -1*b/m
    return vth
        
def vth_trans(i,v):
    v_,gm = deriv(v,i,window=4)
    
    v__,dgm = deriv(v_,gm,window=6)
    
    max_V = v__[np.argmax(dgm)]
    slope = np.max(dgm)
    
    idx = np.where(v_ == max_V)[0][0]
    
    vth = max_V - gm[idx]/slope
    
    return vth

class Model:
    def __init__(self,temperature):
        assert temperature == 300 or temperature == 350, 'Invalid temperature. 300 or 350 only.'
        self.temp = int(temperature)
        self.VTH_LIN_REF = 3.761
        self.VTH_TRANS_REF = 2.671
        if self.temp == 350:
            self.VTH_LIN_REF = 3.248
            self.VTH_TRANS_REF = 2.298
            
        # keep objects for vths incase weighting is desired in the cost function
        self.REF = np.array([self.VTH_LIN_REF,self.VTH_TRANS_REF])
        
    # Sum of squares cost function
    def cost(self,IdVdfile,IdVgfile):
        # input and interpret log file
        # Vth
        Id, Vg, _ = IV_data(IdVgfile)
        vth_lin = vth_linear(Id,Vg)
        vth_tran = vth_trans(Id,Vg)
        # Vbd
        # not yet
        output = np.array([vth_lin,vth_tran])
        
        cost = np.sum(np.power(output-self.REF,2))
        return cost

    # Generates semiconductor geometry
    def buildDeck(self,particle,IdVg,IdVd,particle_num):
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
        # if geometry is bad, flag for ensuring high cost
        if p_source_width <= 0:
            return 1
        else:
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
            
            ATLAS.deck(nsub,ndrift,p,nsource,source,oxide,gate,p_doping,n_drift_doping,n_plus_doping,dit,IdVg,IdVd,filename,boundary,particle_num,tox,self.temp)
            
            return filename


    def simulate(self,particles,run_num):
#        # Testing
#        particle_costs = np.zeros(particles.shape[1])
#        if run_num == 0:
#            # return costs from existing .log files (5)
#            for i,x in enumerate(['Vg_p00.log','Vg_p10.log','Vg_p11.log','Vg_p12.log','Vg_p13.log']):
#                particle_costs[i] = self.cost('','../Data/Model/raw/%s'%x)
#                print '#%d / 5'%(i+1)
#        elif run_num == 1:
#            # return costs from existing .log files (5)
#            for i,x in enumerate(['Vg_p20.log','Vg_p21.log','Vg_p22.log','Vg_p23.log','Vg_p24.log']):
#                particle_costs[i] = self.cost('','../Data/Model/raw/%s'%x)
#                print '#%d / 5'%(i+1)
#        return particle_costs
#        # end testing section
        
        particle_costs = np.zeros(particles.shape[1])
        
        for particle_num in range(particles.shape[1]):
            particle = particles[:,particle_num].flatten()
            IdVdfile = 'SiC_IdVd_run_%d_particle_%d.log'%(run_num,particle_num)
            IdVgfile = 'SiC_IdVg_run_%d_particle_%d.log'%(run_num,particle_num)
            # build deck commands
            deckFile = self.buildDeck(particle,IdVgfile,IdVdfile,particle_num)
            
            if deckFile == 1:
                #penalize arbitrarily high
                particle_costs[particle_num] = 1e15
                print '#%d / %d  -- Invalid geometry'%(particle_num+1,particles.shape[1])
            else:
                # call simulation
                cmd = '\\sedatools\\exe\\deckbuild -run %s'%deckFile
                subprocess.call(cmd.split(' '))
                
                # calculate cost
                particle_costs[particle_num] = self.cost(IdVdfile,IdVgfile)
                print '#%d / %d'%(particle_num+1,particles.shape[1])
        
        return particle_costs

class PSO:
    def __init__(self,inertia,social,cognitive,particleNum,iterationNum,elementNum=6):
        self.NUMBER_OF_PARTICLES = particleNum
        self.NUMBER_OF_ITERATIONS = iterationNum
        self.NUMBER_OF_ELEMENTS = elementNum
        self.INERTIAL_CONST = inertia
        self.SOCIAL_COMP = social
        self.COGNITIVE_COMP = cognitive
        
    def velocityUpdate(self,velocities,particle_vectors, particle_bests_v, global_best_v,r,s):
        # make global_best_V a matrix instead of a vector
        global_best_v = np.repeat(global_best_v,particle_vectors.shape[1],axis=1)
        
        # New velocity with inertial scaling
        vel = np.multiply(velocities,self.INERTIAL_CONST)
        # Add congnitive component to the velocity
        vel = vel + np.multiply(particle_bests_v-particle_vectors,r*self.COGNITIVE_COMP)
        # Add social component to the velocity
        vel = vel + np.multiply(global_best_v-particle_vectors,s*self.SOCIAL_COMP)
        print '...updated velocity'
        return vel

    def optimize(self,Model):
        # CONSTANTS
        
        # particle structure
        #
        # | n+width | Lch | N-drift/gate length | n+ | n- | dit | (in column vector form)
        ##########################################################
        
        print 'Initializing particles...'
        
        particle = np.zeros((self.NUMBER_OF_ELEMENTS,1),dtype=np.float64)
        
        # initialize particles
        particle[0,0] = 182
        particle[1,0] = 195
        particle[2,0] = 200
        particle[3,0] = 5e17
        particle[4,0] = 1e16
        particle[5,0] = 3e10
        #print particle
        
        for x in range(self.NUMBER_OF_PARTICLES):
            if x>0:
                noise = 0.5 + np.random.random(size=(self.NUMBER_OF_ELEMENTS,1))
                particle = np.hstack((particle,np.multiply(particle[:,0].reshape(self.NUMBER_OF_ELEMENTS,1),noise)))
                
        #print particle.shape
        print 'Iteration 0:'
        # initialize velocities
        v = np.random.rand(self.NUMBER_OF_ELEMENTS,self.NUMBER_OF_PARTICLES)*2 - 1
        print 'running simulations...'
        # run initial simulations
        particle_costs = Model.simulate(particle,0)
        
        # save particle-bests (intially equivalent to the first round)
        particle_best_costs = particle_costs.copy()
        particle_best_vectors = particle.copy()
        
        # save global-best
        global_best_cost = np.min(particle_best_costs)
        global_best_vector = particle[:,np.argmin(particle_best_costs)].reshape(self.NUMBER_OF_ELEMENTS,1)
        saveFile = 'logs/global_best_particle_iteration_0.npy'
        np.save(saveFile,global_best_vector.flatten())
        
        costs = [global_best_cost]
        print 'simulations done.\n'
        # main loop
        for q in range(self.NUMBER_OF_ITERATIONS):
            print 'Iteration %d:'%(q+1)
            
            # randomize r and s vectors
            r = np.random.rand(self.NUMBER_OF_ELEMENTS,self.NUMBER_OF_PARTICLES)
            s = np.random.rand(self.NUMBER_OF_ELEMENTS,self.NUMBER_OF_PARTICLES)
            
            # update velocities
            v = self.velocityUpdate(v, particle, particle_best_vectors, global_best_vector,r,s)
            
            # update particles
            particle = particle + v
            print 'particles updated...'
            print 'running simulations...'
            # run simulation for each particle
            particle_costs = Model.simulate(particle,q+1)
            
            # save particle bests
            for i,particle_cost in enumerate(particle_costs):
                if particle_cost < particle_best_costs[i]:
                    particle_best_costs[i] = particle_cost
                    particle_best_vectors[:,i] = particle[:,i]
            
            # save global best
            global_best_cost = np.min(particle_best_costs)
            costs.append(global_best_cost)
            global_best_vector = particle_best_vectors[:,np.argmin(particle_best_costs)].reshape(self.NUMBER_OF_ELEMENTS,1)
            saveFile = 'logs/global_best_particle_iteration_%d.npy'%(q+1)
            np.save(saveFile,global_best_vector.flatten())
            
            #save particle positions
            saveFile = 'logs/particle_iteration_%d.npy'%(q+1)
            np.save(saveFile, particle)
            
            print 'simulations done.\n'
        
        #save global best cost history
        np.save('logs/global_best_costs.npy',np.array(costs))
        print 'Done...'
        print 'Best cost: %.2f'%costs[-1]
        print 'Best Particle:'
        for x in global_best_vector.flatten():
            print x

def main(argv):
    # create and maintain logs
    import os, time
    if 'logs' not in os.listdir('.'):
        os.mkdir('logs')
        time = '_'.join([str(x) for x in time.gmtime()[:-4]])
        os.system('echo %s >> ./logs/time.txt'%time)
    else:
        with open('./logs/time.txt','rb') as fObj:
            ti = fObj.read()
        os.system('mv logs old_logs_%s'%ti)
        os.system('mkdir logs')
        time = '_'.join([str(x) for x in time.gmtime()[:-4]])
        os.system('echo %s >> ./logs/time.txt'%time)
    
    number_of_particles = 0
    number_of_iterations = 0
    ambient_temperature = 300
    inertial_component = 0
    social_component = 0
    cognitive_component = 0
    help_flag = False
    while argv:
        if argv[0] == '--inertial':
            inertial_component = float(argv[1])
        elif argv[0] == '--social':
            social_component = float(argv[1])
        elif argv[0] == '--cognitive':
            cognitive_component = float(argv[1])
        elif argv[0] == '--particles':
            number_of_particles = int(argv[1])
        elif argv[0] == '--iterations':
            number_of_iterations = int(argv[1])
        elif argv[0] == '--temp':
            ambient_temperature = int(argv[1])
        elif argv[0] == '-h' or argv[0] == '--help':
            help_flag = True
        argv = argv[1:]
    
    if help_flag:
        print ''
        print ' --inertial <inertial component>'
        print ' --social <social component>'
        print ' --cognitive <cognitive component>'
        print ' --particles <particle number>'
        print ' --iterations <iteration number>'
        print ' --temp <300 or 350> (optional, 300 by default)'
    elif not (number_of_iterations and number_of_particles and inertial_component and social_component and cognitive_component):
        print 'Flag missing!'
        print ' --inertial <inertial component>'
        print ' --social <social component>'
        print ' --cognitive <cognitive component>'
        print ' --particles <particle number>'
        print ' --iterations <iteration number>'
        print ' --temp <300 or 350> (optional, 300 by default)'
    else:
        MosfetModel = Model(ambient_temperature)
        Swarm = PSO(inertial_component, social_component, cognitive_component, number_of_particles, number_of_iterations)
        Swarm.optimize(MosfetModel)

if __name__ == '__main__':
    from sys import argv
    main(argv)