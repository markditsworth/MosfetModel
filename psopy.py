#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 17:30:16 2018

@author: markditsworth
"""

import numpy as np
import matplotlib.pyplot as plt

class PSO:
    def __init__(self, user_fxn, initial, swarm_size=20, stopping_criteria='iteration', distance_tolerance=None,
                 error_tolerance=0.01, iter_number=100, inertia=0.9, cognitive=2.0, social=2.0, method='standard',
                 logging=False, max_iter=np.inf, swarm_build=None, report=False, visual_dim=(0,1)):
        # ensure kwargs are valid
        assert swarm_build in ['+-50%','random',None], "Invalid swarm_build.\n'+-50%', 'random', or None" 
        assert stopping_criteria.lower() in ['error','iteration','distance'], "Invalid stopping criteria.\n 'iteration' (default),'error', or 'distance'"
        assert method in ['standard'], "Invalide method." #add more later
        # set basic PSO object values
        self.cost_fcn = user_fxn
        self.inertia = inertia
        self.cognitive = cognitive
        self.social = social
        self.initial = initial
        
        # handle stopping criteria
        self.stopping = stopping_criteria.lower()
        self.tol = iter_number
        if self.stopping == 'error':
            self.tol = error_tolerance
        elif self.stopping == 'distance':
            self.tol = distance_tolerance
        assert not self.tol == None, "tolerance cannot be NoneType. distance_tolerance must be set."
        
        self.max_iter = max_iter
        if self.max_iter != np.inf and self.stopping == 'iteration':
            print 'Stopping criteria set to be %d iterations. Ignoring max_iter setting.'%self.tol
            self.max_iter = np.inf
        
        # PSO method, and swarm building settings
        self.method = method
        self.swarm_size = swarm_size
        self.swarm_build = swarm_build
        # if initial swarm particles are provided, use the number provided as swarm size
        # otherwise, if an array is given, a swarm must be built using the swarm_build parameter
        if len(self.initial.shape) > 1:
            self.swarm_size = self.initial.shape[1]
            if swarm_build != None:
                print 'Swarm provided. swarm_build will be ignored.'
                self.swarm_build = None
                
        self.logging = logging
        self.report = report
        self.visual_dim = visual_dim
        
    
    def checkStop(self,iteration,particles,cost):
        if self.stopping == 'error':
            if cost <= self.tol:
                return True
        elif self.stopping == 'distance':
            print 'distance not accounted for yet'
            return True
        else:
            if iteration > self.tol:
                return True
        
        return False
    
    def makeParticles(self):
        # generate additional particles with value varying between 50% more or less than initial values
        if self.swarm_build == '+-50%':
            # make column vector
            particles = np.reshape(self.initial,(np.shape(self.initial)[0],1))
            # add columns to build to desired swarm size
            for x in range(self.swarm_size):
                if x>0:
                    #randomly vary between 50% and 150% of the initial values
                    noise = 0.5 + np.random.random(size=(particles.shape[0],1))
                    particles = np.hstack((particles,np.multiply(particles[:,0].reshape(self.NUMBER_OF_ELEMENTS,1),noise)))
            return particles
        # generate additional particles vith values varying randomly between the minimum and maximum initial values
        elif self.swarm_build == 'random':
            # make column vector
            particles = np.reshape(self.initial,(np.shape(self.initial)[0],1))
            # add columns to build to desired swarm size
            for x in range(self.swarm_size):
                if x>0:
                    #randomly vary between min and max of the initial values
                    diff = np.max(self.initial) - np.min(self.initial)
                    new = (np.random.random(size=(particles.shape[0],1)) * diff) + np.min(self.initial)
                    particles = np.hstack((particles,new))
            return particles
        # the initial input is the swarm
        else:
            return self.initial.copy()
    
    def standardVelocityUpdate(self,particles,velocities,pbs,g):
        #ensure g is a column vector
        g = g.reshape((particles.shape[0],1))
        
        # create random vectors
        r = np.random.rand(particles.shape[0],particles.shape[1])
        s = np.random.rand(particles.shape[0],particles.shape[1])
        
        # New velocity with inertial scaling
        vel = np.multiply(velocities,self.inertia)
        
        # Add congnitive component to the velocity
        vel = vel + np.multiply(pbs-particles, r*self.cognitive)
        
        # Add social component to the velocity
        vel = vel + np.multiply(g.repeat(particles.shape[1],axis=1)-particles, s*self.social)
        
        return vel
    
    def swarm(self,verbose=False):
        #initialize particles
        particles = self.makeParticles()
        
        #initialize random velocities [-1,1]
        velocities = np.random.rand(particles.shape[0],particles.shape[1])*2 - 1
        
        # initialize pbs and g arbitrarily high
        particle_best_costs = np.ones(self.swarm_size) * 1e20
        global_best_cost = 1e20
        
        # initialize best_particle and particle_bests
        particle_bests = particles.copy()
        best_particle = np.zeros(particles.shape[0])
        
        #housekeeping
        iteration_number = 0
        
        cost_history = np.array([])
        while True:
                
            # Break from loop if stopping criteria met
            if self.checkStop(iteration_number,particles,global_best_cost):
                print 'Stopping criteria reached! %d iterations.'%(iteration_number-1)
                break
            if iteration_number >= self.max_iter:
                print 'WARNING: No convergence in %d iterations. Results not valid.'%self.max_iter
                break
            
            if verbose: print '\niteration %d....'%iteration_number
            
            if self.logging:
                self.save_swarm(particles,iteration_number)
            
            # Run calculation
            costs = self.cost_fcn(particles)
            
            if verbose: print 'costs calculated...'
            
            #update particle_best_costs
            for i,cost in enumerate(costs):
                if cost < particle_best_costs[i]:
                    particle_best_costs[i] = cost
                    
                    #update particle_bests
                    particle_bests[:,i] = particles[:,i]
            
            if verbose:
                print '\nparticle bests updated:'
                for i,x in enumerate(particle_best_costs):
                    print 'particle %d cost: %f'%(i,x)
            
            #update global_best_cost
            global_best_cost = min(global_best_cost,np.min(costs))
            cost_history = np.hstack((cost_history,global_best_cost))
            
            if verbose: print '\nglobal best cost: %f'%global_best_cost
            
            #update best_particle
            best_particle = particles[:,np.argmin(costs)]
            
            if verbose:
                print '\nbest particle found:'
                for x in best_particle:
                    print x
            
            #update velocities
            if self.method == 'standard':
                velocities = self.standardVelocityUpdate(particles,velocities,particle_bests,best_particle)
            else:
                print 'No other methods are available at this time.'
                break
            
            #update particles
            particles = np.add(particles,velocities)
            iteration_number += 1
            
        if self.logging:
            np.save('cost_hist.npy',cost_history)
            if self.report:
                self.graph(self.visual_dim,iteration_number-1)
            
        return best_particle, global_best_cost
    
    def save_swarm(self,swarm,iteration_count):
        fname = 'particles_iteration_%d.npy'%iteration_count
        np.save(fname,swarm)
        
    
    def graph(self,pair,iterations):
        # handle large number of iterations, only show 20
        steps = np.arange(iterations)
        if iterations>20:
            step_size = iterations/20
            offset = iterations % 20
            steps = np.arange(offset, 20*step_size, step_size)
            
        # 1x2 figure
        fig,axes = plt.subplots(2)
        # for each of the selected 20 iterations
        for i in steps:
            fname = 'particles_iteration_%d.npy'%i
            swarm = np.load(fname)
            ### SELECT COLOR FOR THIS ITERATION
            c=np.random.random(3)
            # plot each particle position for that iteration
            flag = 1
            for j in range(swarm.shape[1]):
                particle = swarm[:,j]
                x = particle[pair[0]]
                y = particle[pair[1]]
                if flag:
                    axes[0].scatter(x,y,s=1,color=c,label='iter. %d'%i) #add color and label later
                else:
                    axes[0].scatter(x,y,s=1,color=c)
                flag = 0
        # axis info
        axes[0].set_title('Particle Swarm')
        #axes[0].legend()
        
        #plot costs
        costs = np.load('cost_hist.npy')
        axes[1].semilogy(costs,c='k')
        # axis info
        axes[1].set(xlabel='iteration number',ylabel='global best cost')
        
        plt.show()
            
            