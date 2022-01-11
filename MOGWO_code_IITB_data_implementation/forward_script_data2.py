
#creating subsurface model
from srfpython import *

ff=open('forward_data.txt','r')
xx=[]
for x in ff:
    xx.append(float(x))
ff.close()

ztop_for=xx[0:5]
vp_for=xx[5:10]
vs_for=xx[10:15]
rh_for=xx[15:20]

dm_for = depthmodel_from_arrays(ztop_for, vp_for, vs_for, rh_for) 
dm_for.write96('iitb_forward_first.mod96')
    
