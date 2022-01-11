
ff=open('iitb_fundmode_all4tech.txt','r')
period_time=[]     
vel_phase=[]
for xxx in ff:
    period_time.append(float(xxx.split()[5]))
    vel_phase.append(float(xxx.split()[6]))
ff.close()
    
fo=open('observed_fund_phase_velocity.txt','w');
for cvel in vel_phase:
    fo.write("%f\n"%(cvel))
fo.close()
