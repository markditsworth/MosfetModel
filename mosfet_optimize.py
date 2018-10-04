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
from psopy import Swarm

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
        self.run_num = 0
        
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
    def buildDeck(self,particle,IdVg,IdVd,particle_num,n_height=0.15,p_height=0.3,total_height=6):
        n_width = particle[0]
        p_width = particle[1]
        n_drift_width = particle[2]
        n_plus_doping = particle[3]
        n_drift_doping = particle[4]
        dit = particle[5]
        
        TOTAL_WIDTH = 600
        TOTAL_HEIGHT = total_height
        DEPTH = 2000
        
        NSUB_HEIGHT = 3
        
        P_HEIGHT = p_height
        
        N_HEIGHT = n_height
        
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
            nh=nHeight,ph=pHeight
            ndrift = "%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f"%(0,TOTAL_HEIGHT-NSUB_HEIGHT, TOTAL_WIDTH,TOTAL_HEIGHT-NSUB_HEIGHT, TOTAL_WIDTH,0, TOTAL_WIDTH-n_drift_width,0, TOTAL_WIDTH-n_drift_width,P_HEIGHT, 0,P_HEIGHT)
            
            p = "%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f"%(0,P_HEIGHT, TOTAL_WIDTH-n_drift_width,P_HEIGHT, TOTAL_WIDTH-n_drift_width,0, TOTAL_WIDTH-n_drift_width-p_width,0, TOTAL_WIDTH-n_drift_width-p_width,N_HEIGHT, p_source_width,N_HEIGHT, p_source_width,0, 0,0)
            
            nsource = "%f,%f %f,%f %f,%f %f,%f %f,%f"%(p_source_width,N_HEIGHT, p_source_width+n_width,N_HEIGHT, p_source_width+n_width,0, source_contact_width,0, p_source_width,0)
            
            source = "%f,%f %f,%f %f,%f %f,%f %f,%f"%(0,0, source_contact_width,0, source_contact_width,-tox, source_contact_width,-0.3, 0,-0.3)
            
            oxide = "%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f"%(source_contact_width,0, p_source_width+n_width,0, TOTAL_WIDTH-n_drift_width,0, TOTAL_WIDTH,0, TOTAL_WIDTH,-tox, p_source_width+n_width-1,-tox, source_contact_width,-tox)
            
            gate = "%f,%f %f,%f %f,%f %f,%f"%(p_source_width+n_width-1,-tox, TOTAL_WIDTH,-tox, TOTAL_WIDTH,-tox-0.3, p_source_width+n_width-1,-tox-0.3)
            
            filename = "SiC_particle_%d.in"%particle_num
            
            ATLAS.deck(nsub,ndrift,p,nsource,source,oxide,gate,p_doping,n_drift_doping,n_plus_doping,dit,IdVg,IdVd,filename,boundary,particle_num,tox,self.temp)
            
            return filename
    
    def _iterate(self):
        n_w = 182
        p_w = 187
        nd_w = 200
        ns = 5e17
        nd = 1e16
        dit=3e10
        p = np.array([n_w,p_w,nd_w,ns,nd,dit])
        num=1
        for TH in [6,7,8]:
            idvd = 'SiC_IdVd_curvatureTest_%d.log'%num
            idvg = 'SiC_IdVg_curvatureTest_%d.log'%num
            deckFile = self.buildDeck(p,idvg,idvd,num,total_height=TH)
            print 'height: %d'%TH
            cmd = '\\sedatools\\exe\\deckbuild -run %s'%deckFile
            subprocess.call(cmd.split(' '))
            num += 1

    def simulate(self,particles):
#        # Testing
#        particle_costs = np.zeros(particles.shape[1])
#        if self.run_num == 0:
#            # return costs from existing .log files (5)
#            for i,x in enumerate(['Vg_p00.log','Vg_p10.log','Vg_p11.log','Vg_p12.log','Vg_p13.log']):
#                particle_costs[i] = self.cost('','../Data/Model/iterations/Iteration2_doping/%s'%x)
#                print '#%d / 5'%(i+1)
#        elif self.run_num == 1:
#            # return costs from existing .log files (5)
#            for i,x in enumerate(['Vg_p20.log','Vg_p21.log','Vg_p22.log','Vg_p23.log','Vg_p24.log']):
#                particle_costs[i] = self.cost('','../Data/Model/iterations/Iteration2_doping/%s'%x)
#                print '#%d / 5'%(i+1)
#        self.run_num += 1
#        return particle_costs
#        # end testing section
        
        particle_costs = np.zeros(particles.shape[1])
        
        for particle_num in range(particles.shape[1]):
            particle = particles[:,particle_num].flatten()
            IdVdfile = 'SiC_IdVd_run_%d_particle_%d.log'%(self.run_num,particle_num)
            IdVgfile = 'SiC_IdVg_run_%d_particle_%d.log'%(self.run_num,particle_num)
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
        
        self.run_num += 1
        
        return particle_costs


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
        #handle OS differences for mv/move
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
    test_flag = False
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
        elif argv[0] == '--test':
            test_flag = True
        argv = argv[1:]
    
    if help_flag:
        print ''
        print ' --inertial <inertial component>'
        print ' --social <social component>'
        print ' --cognitive <cognitive component>'
        print ' --particles <particle number>'
        print ' --iterations <iteration number>'
        print ' --temp <300 or 350> (optional, 300 by default)'
        print ' --test Runs the _iterate function to test conditions'
    elif not (number_of_iterations and number_of_particles and inertial_component and social_component and cognitive_component):
        print 'Flag missing!'
        print ' --inertial <inertial component>'
        print ' --social <social component>'
        print ' --cognitive <cognitive component>'
        print ' --particles <particle number>'
        print ' --iterations <iteration number>'
        print ' --temp <300 or 350> (optional, 300 by default)'
    elif test_flag:
        MosfetModel = Model(ambient_temperature)
        MosfetModel._iterate()
    else:
        MosfetModel = Model(ambient_temperature)
        # initialize swarm
        particle = np.zeros((6,number_of_particles),dtype=np.float64)
        particle[0,0] = 182
        particle[1,0] = 195
        particle[2,0] = 200
        particle[3,0] = 5e17
        particle[4,0] = 1e16
        particle[5,0] = 3e10
        #print particle
        
        for x in range(number_of_particles):
            if x>0:
                particle[0,x] = np.random.randint(160,200)
                particle[1,x] = np.random.randint(160,200)
                particle[2,x] = np.random.randint(160,200)
                particle[3,x] = np.random.randint(1,100) *np.power(10,16)
                particle[4,x] = np.random.randint(1,100) *np.power(10,15)
                particle[5,x] = np.random.randint(1,100) *np.power(10,9)
                
                #noise = 0.5 + np.random.random(size=(6,1))
                #particle = np.hstack((particle,np.multiply(particle[:,0].reshape(6,1),noise)))
        
        #pso = PSO(MosfetModel.simulate, particle, inertia=inertial_component, social=social_component,
        #          cognitive=cognitive_component, iter_number=30, logging=True)
        #best,cost = pso.swarm(verbose=True)
        swarm = Swarm(MosfetModel.simulate, particle, stopping_criteria='iteration', iter_number=30,
                      logging=True)
        # gbest
        best,cost = swarm.gbest(inertial=inertial_component, social=social_component, cognitive=cognitive_component,
                                verbose=True)
        # lbest
        '''
        best,cost = swarm.lbest(inertial=inertial_component, social=social_component, cognitive=cognitive_component,
                                n_neighbors=3,verbose=True)
        # FIPS
        best,cost = swarm.FIPS()
        
        # FDR
        best,cost = swarm.FDR(inertial=inertial_component, social=social_component, cognitive=cognitive_component,
                              verbose=True)'''
        
        print 'best cost: %f'%cost
        print best

if __name__ == '__main__':
    from sys import argv
    main(argv)