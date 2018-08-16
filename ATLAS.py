#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 12:57:12 2018

@author: markditsworth
"""

def deck(nsub,ndrift,p,n_plus,source,oxide,gate,p_doping,n_drift_doping,n_plus_doping,dit,IdVg_fname,IdVd_fname,filename,boundary,num,tox):
    deck = r"""
# (c) Silvaco Inc., 2015
go devedit

DevEdit version=2.4.0.R

work.area x1=0 y1=-0.68 x2=10 y2=15
# devedit 2.4.0.R (Thu May  8 12:10:27 PDT 1997)
# libsflm 2.0.0.R (Thu May  1 18:03:38 PDT 1997)
# libDW_Misc 1.20.0.R (Mon Apr 28 17:55:25 PDT 1997)
# libCardDeck 1.20.0.R (Tue Apr 29 15:01:54 PDT 1997)
# libGeometry 1.20.0.R (Mon Apr 28 18:17:55 PDT 1997)
# libDW_Set 1.20.0.R (Mon Apr 28 17:57:52 PDT 1997)
# libSVC_Misc 1.20.0.R (Mon Apr 28 18:20:53 PDT 1997)
# libSDB 1.0.6.C (Mon May  5 16:28:49 PDT 1997)
# libSSS 1.20.0.R (Mon May  5 16:29:45 PDT 1997)
# libMeshBuild 1.20.0.R (Wed May  7 23:57:48 PDT 1997)
# libDW_Make 1.1.3.R (Thu May  1 20:07:31 PDT 1997)

region reg=1 name=n+sub mat=4H-SiC color=0x7f00ff pattern=0x8 \
	polygon="%s"
#
impurity id=1 region.id=1 imp=Arsenic \
	peak.value=1e19 ref.value=1000000000000 comb.func=Multiply
#
constr.mesh region=1 default max.height=2 min.width=100

region reg=2 name=n-drift mat=4H-SiC color=0x7f00ff pattern=0x8 \
	polygon="%s"
#
impurity id=1 region.id=2 imp=Arsenic \
	peak.value=%f ref.value=1000000000000 comb.func=Multiply
#
constr.mesh region=2 default max.height=1.25 min.width=100

region reg=3 name=pbase mat=4H-SiC color=0x7f00ff pattern=0x8 \
	polygon="%s"
#
impurity id=1 region.id=3 imp=Boron \
	peak.value=%f ref.value=1000000000000 comb.func=Multiply
#
constr.mesh region=3 default min.width=100

region reg=4 name=n+source mat=4H-SiC color=0x7f00ff pattern=0x8 \
	polygon="%s"
#
impurity id=1 region.id=4 imp=Arsenic \
	peak.value=%f ref.value=1000000000000 comb.func=Multiply
#
constr.mesh region=4 default min.width=100
#max.width=0.25

region reg=5 name=ox mat="Silicon Oxide" \
	polygon="%s"
#
constr.mesh region=5 default min.width=100 
#max.height=0.05

region reg=6 name=source mat=Aluminum elec.id=1 work.func=0 color=0xffc8c8 pattern=0x7 \
	polygon="%s"
#
constr.mesh region=6 default min.width=100 
#max.height=0.1

region reg=7 name=gate mat=Aluminum elec.id=2 work.func=0 color=0xffc8c8 pattern=0x7 \
	polygon="%s"
#
constr.mesh region=7 default min.width=100 
#max.height=0.1

substrate name="drain" electrode=3 workfunction=0


# Set Meshing Parameters
#
base.mesh height=5 width=5
#
bound.cond !apply max.slope=30 max.ratio=100 rnd.unit=0.001 line.straightening=1 align.points when=automatic
#
imp.refine imp="Net Doping" scale=log transition=1e+10
imp.refine min.spacing=0.02
#
#constr.mesh  max.angle=90 max.ratio=300 max.height=1000 \
#	max.width=1000 min.height=0.0001 min.width=0.0001
constr.mesh  max.angle=90 max.ratio=300 max.height=1000 \
	max.width=1000 min.height=0.0001 min.width=4
#
constr.mesh type=Semiconductor default 
#min.width=4
#
constr.mesh type=Insulator default
#
constr.mesh type=Metal default
#
constr.mesh type=Other default
#
constr.mesh region=1 default
#
constr.mesh region=2 default
#
constr.mesh region=3 default 
#min.width=4
#
#constr.mesh region=4 default max.width=0.25
constr.mesh region=4 default 
#min.width=4

#SiO2
constr.mesh region=5 default max.height=%f max.width=4
#min.width=4 max.height=0.05

#Source
#constr.mesh region=6 default max.width=0.25
constr.mesh region=6 default max.height=0.05 min.width=4
#min.width=5 max.height=0.1

#Gate
constr.mesh region=7 default max.height=0.1 max.width=4
#max.height=0.1
##########################################
constr.mesh id=1 under.reg=ox depth=0.1 default max.height=0.0125 max.width=0.4
##constr.mesh id=1 under.reg=gate depth=0.1 default max.height=0.03 max.width=5

##constr.mesh id=3 under.reg=source depth=0.1 default max.height=0.02 max.width=0.25
constr.mesh id=3 under.reg=source depth=0.1 default max.height=0.0125 max.width=4

#HERE
#constr.mesh id=2 x1=0 x2=500 y1=0 y2=2   default max.height=0.2  min.width=4
#constr.mesh id=6 x1=0 x2=400 y1=0 y2=0.8 default max.height=0.1  min.width=4
#constr.mesh id=1 x1=0 x2=350 y1=0 y2=0.1 default max.height=0.05 max.width=2

#constr.mesh id=5 x1=0 x2=600 y1=4 y2=11 default min.width=25 min.height=2
# N/N-drift boundary
constr.mesh id=4 x1=0 y1=%s x2=1e+06 y2=%s default max.height=0.2 min.width=4

##constr.mesh id=5 x1=0 y1=2 x2=1e+06 y2=4 default max.height=0.6
Mesh Mode=MeshBuild


base.mesh height=5 width=5

bound.cond !apply max.slope=30 max.ratio=100 rnd.unit=0.001 line.straightening=1 align.Points when=automatic

struct outf=Wolfspeed_SiC_%d.str
#tonyplot Wolfspeed_SiC.str
##############################

#go atlas
#
#mesh infile=Wolfspeed_SiC_%d.str
#
#material material=4H-SiC permitti=9.66 eg300=3.26 \
#			edb=0.1 gcb=2 eab=0.2 gvb=4 \
#			nsrhn=3e17 nsrhp=3e17 taun0=5e-10 taup0=1e-10 \
#			tc.a=100 taurel.el=2.8e-12
#models temp=300 conmob bgn srh print
#impact selb
#mobility material=4H-SiC	vsatn=2e7 vsatp=2e7 betan=2 betap=2 \
#			mu1n.caug=10  mu2n.caug=410 ncritn.caug=13e17  \
#			deltan.caug=0.6 gamman.caug=0.0 \
#			alphan.caug=-3 betan.caug=-3 \
#			mu1p.caug=20   mu2p.caug=95  ncritp.caug=1e19 \
#			deltap.caug=0.5  gammap.caug=0.0 \
#			alphap.caug=-3 betap.caug=-3
#
#contact name=gate n.poly
##contact name=source workfunction=0
##contact name=drain workfunction=0
#
##interface qf=3e10
#interface qf=%f
##
#solve init
#method newton trap maxtraps=20 climit=1e-4  ir.tol=1e-12 ix.tol=1e-12
#
# 
#solve init
##
#log outf=%s
##solve vgate=0
#
#solve vsource=0
#solve vdrain=0.03
#solve vdrain=0.05
#solve vdrain=0.1
#solve vgate=0
#solve vdrain=0.09
#solve vdrain=0.25 vstep=0.25 vfinal=5 name=drain 
#solve vstep=1.0 vfinal=10 name=drain
#solve vstep=10.0 vfinal=1000 name=drain
#solve vstep=100.0 vfinal=3000 name=drain

##############################

go atlas

mesh infile=Wolfspeed_SiC_%d.str

material material=4H-SiC permitti=9.66 eg300=3.26 \
			edb=0.1 gcb=2 eab=0.2 gvb=4 \
			nsrhn=3e17 nsrhp=3e17 taun0=5e-10 taup0=1e-10 \
			tc.a=100 taurel.el=2.8e-12
models temp=300 conmob bgn srh print
#impact selb
mobility material=4H-SiC	vsatn=2e7 vsatp=2e7 betan=2 betap=2 \
			mu1n.caug=10  mu2n.caug=410 ncritn.caug=13e17  \
			deltan.caug=0.6 gamman.caug=0.0 \
			alphan.caug=-3 betan.caug=-3 \
			mu1p.caug=20   mu2p.caug=95  ncritp.caug=1e19 \
			deltap.caug=0.5  gammap.caug=0.0 \
			alphap.caug=-3 betap.caug=-3

contact name=gate n.poly
#contact name=source workfunction=0
#contact name=drain workfunction=0

#interface qf=3e10
interface qf=%f
#
solve init
method newton trap maxtraps=20 climit=1e-4  ir.tol=1e-12 ix.tol=1e-12

 
solve init
#
log outf=%s
#solve vgate=0

solve vsource=0
solve vdrain=0.03
solve vdrain=0.05
solve vdrain=0.1
solve vgate=0.01
solve vgate=0.05 vstep=0.05 vfinal=8 name=gate

 
#tonyplot
quit
    """%(nsub,ndrift,n_drift_doping,p,p_doping,n_plus,n_plus_doping,oxide,source,gate,tox/2.0,boundary-1,boundary+1,num,num,dit,IdVd_fname,num,dit,IdVg_fname)
    
    # save deck as
    with open(filename,'wb') as fObj:
        fObj.write(deck)
    
    
    