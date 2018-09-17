#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 17:30:16 2018

@author: markditsworth
"""

import numpy as np
import matplotlib.pyplot as plt

class Swarm:
    def __init__(self, user_fxn, initial, swarm_size=20, stopping_criteria='iteration', distance_tolerance=None,
                 error_tolerance=0.01, iter_number=100,logging=False, max_iter=np.inf, swarm_build=None,
                 report=False, visual_dim=(0,1)):
        # ensure kwargs are valid
        assert swarm_build in ['+-50%','random',None], "Invalid swarm_build.\n'+-50%', 'random', or None" 
        assert stopping_criteria.lower() in ['error','iteration','distance'], "Invalid stopping criteria.\n 'iteration' (default),'error', or 'distance'"
        
        # set basic PSO object values
        self.cost_fcn = user_fxn
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
    
    def standardVelocityUpdate(self,particles,velocities,pbs,g,inertial,social,cog):
        #ensure g is a column vector
        g = g.reshape((particles.shape[0],1))
        
        # create random vectors
        r = np.random.rand(particles.shape[0],particles.shape[1])
        s = np.random.rand(particles.shape[0],particles.shape[1])
        
        # New velocity with inertial scaling
        vel = np.multiply(velocities,inertial)
        
        # Add congnitive component to the velocity
        vel = vel + np.multiply(pbs-particles, r*cog)
        
        # Add social component to the velocity
        vel = vel + np.multiply(g.repeat(particles.shape[1],axis=1)-particles, s*social)
        
        return vel
    
    def lbestVelocityUpdate(self,particles,velocities,pbs,pbest_costs,kdtree,inertial,social,cognitive,K):
        # create random vectors
        r = np.random.rand(particles.shape[0],particles.shape[1])
        s = np.random.rand(particles.shape[0],particles.shape[1])
        
        # New velocity with inertial scaling
        vel = np.multiply(velocities,inertial)
        
        # Add congnitive component to the velocity
        vel = vel + np.multiply(pbs-particles, r*cognitive)
        
        # Add social component to the velocity from K nearest neighbors
        for p in range(particles.shape[1]):
            _, nn_idx = kdtree.query(particles[:,p],k=K+1) # get indeces of K nearest neighbors (but include itself)
            neighborhood_best_idx = nn_idx[np.argmin(pbest_costs[nn_idx])] # get index of neighbor with lowest cost
            vel[:,p] = vel[:,p] + np.multiply(pbs[:,neighborhood_best_idx] - particles[:,p],s[:,p]*social)
        
        return vel
    
    def fipsVelocityUpdate(self,particles,velocities,pbs,phi,chi,A):
        vel = np.zeros(particles.shape)
        if type(A) == str:
            if A == 'ones':
                A = np.ones((particles.shape[1],particles.shape[1]))
        
        # create random vector
        r = phi*np.random.rand(particles.shape[0],particles.shape[1])
        
        for i in range(particles.shape[1]):
            diffs = pbs - particles[:,i].reshape((particles.shape[0],1))
            vel[:,i] = np.sum(np.multiply(A[:,i],diffs), axis=1)
        
        vel = ((np.multiply(r,vel) / particles.shape[1]) + velocities)*chi
        
        return vel
        
    
    def optimize(self,verbose,method,kwargs={}):
        if method == 'lbest':
            from scipy.spatial import KDTree
            
        elif method == 'binary':
            # input should only be ones and zeros
            check = (self.initial == np.ones(self.initial.shape)).astype(int) + (self.initial == np.zeros(self.initial.shape)).astype(int)
            assert check.all() == True, 'Invalid initial swarm. For binary PSO, particles must be in [0,1]'
            
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
            
            if self.logging:
                self.save_swarm(particles,iteration_number)
                
            if method == 'lbest':
                kdt = KDTree(particles.T) #KD-Tree built on row vectors, not column
            
            # Run calculation
            costs = self.cost_fcn(particles)
            
            #update particle_best_costs
            for i,cost in enumerate(costs):
                if cost < particle_best_costs[i]:
                    particle_best_costs[i] = cost
                    
                    #update particle_bests
                    particle_bests[:,i] = particles[:,i]
            
            #update global_best_cost
            global_best_cost = min(global_best_cost,np.min(costs))
            cost_history = np.hstack((cost_history,global_best_cost))
            
            
            #update best_particle
            best_particle = particles[:,np.argmin(costs)]
            
            if verbose:
                print 'iteration %d...'%iteration_number
                for i,x in enumerate(particle_best_costs):
                    if x == global_best_cost:
                        print 'particle %d cost: %.3e <-- GLOBAL BEST'%(i,x)
                    else:
                        print 'particle %d cost: %.3e'%(i,x)
            else:
                print 'iteration %6d | best cost = %.3e'%(iteration_number,global_best_cost)
            
            #update velocities
            if method == 'standard' or method == 'binary':
                velocities = self.standardVelocityUpdate(particles,velocities,particle_bests,best_particle,
                                                         kwargs['inertial'],kwargs['social'],kwargs['cognitive'])
            elif method == 'FIPS':
                velocities = self.fipsVelocityUpdate(particles,velocities,particle_bests,
                                                     kwargs['phi'],kwargs['chi'],kwargs['adjacency'])
            elif method == 'lbest':
                velocities = self.lbestVelocityUpdate(particles,velocities,particle_bests,particle_best_costs,
                                                      kdt,kwargs['inertial'],kwargs['social'],kwargs['cognitive'],
                                                      kwargs['K'])
            else:
                print 'No other methods are available at this time.'
                break
            
            #update particles
            if method == 'binary':
                probabilities = 1/(1 + np.exp(-1*velocities)) # velocities to [0,1] via sigmoid function
                r = np.random.random(particles.shape)
                particles = (r < probabilities).astype(int)
            else:
                particles = np.add(particles,velocities)
                
            iteration_number += 1
            
        if self.logging:
            np.save('cost_hist.npy',cost_history)
            if self.report:
                self.graph(self.visual_dim,iteration_number-1)
            
        return best_particle, global_best_cost
    
    def gbest(self,inertial=0.6,social=1.7,cognitive=0.4,verbose=False):
        kwargs = {'inertial':inertial,'social':social,'cognitive':cognitive}
        return self.optimize(verbose,'standard',kwargs=kwargs)
    
    def FIPS(self,A='ones',phi=4.1,chi=0.7,verbose=False):
        kwargs = {'phi':phi,'chi':chi,'adjacency':A}
        return self.optimize(verbose,'FIPS',kwargs=kwargs)
    
    def binary(self,inertial=0.6,social=1.7,cognitive=0.4,verbose=False):
        kwargs = {'inertial':inertial,'social':social,'cognitive':cognitive}
        return self.optimize(verbose,'binary',kwargs=kwargs)
    
    def lbest(self,inertial=0.6,social=1.7,cognitive=0.4,n_neighbors=2,verbose=False):
        kwargs = {'inertial':inertial,'social':social,'cognitive':cognitive,'K':n_neighbors,
                  'verbose':verbose}
        return self.optimize(verbose,'lbest',kwargs=kwargs)
    
    def FDR(self):
        return 'Not yet implemented'
    
    def save_swarm(self,swarm,iteration_count):
        fname = 'particles_iteration_%d.npy'%iteration_count
        np.save(fname,swarm)
        
    def graph(self,pair,iterations):
        # handle large number of iterations, only ever show 20
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
            
            