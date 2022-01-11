
#Storing first mode phase velocities
ff=open('iitb_forward_first.surf96','r')
period_time_for=[]      
vel_phase_for=[]
for xxx in ff:
    period_time_for.append(float(xxx.split()[5]))
    vel_phase_for.append(float(xxx.split()[6]))
ff.close()
    
fo=open('forward_first_phase_velocity.txt','w');
for cvel in vel_phase_for:
    fo.write("%f\n"%(cvel))
fo.close()
