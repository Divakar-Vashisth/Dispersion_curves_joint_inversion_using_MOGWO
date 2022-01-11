function [out]= joint_inv_obj_fun(x,vel_phase_fund,vel_phase_first)

d=[0.0,x(1),x(2),x(3),x(4)]./1000; %conversion to Km
vp=x(5:9)'./1000; %conversion to Km/s
vs=x(10:14)'./1000; %conversion to Km/s
ro=x(15:19)'; %if not inverting for density, replace with ro=[ro1,ro2...] where ro1, ro2 are assumed or known densities of first and second layer respectively


final_x=[d,vp,vs,ro];


dlmwrite('forward_data.txt',final_x,'\n')

[status,cmdout]=unix ("forward_script1.scr"); %for fundamental mode
vel_phase_fund_for=load('forward_fund_phase_velocity.txt');

if length(vel_phase_fund_for)~= length(vel_phase_fund)
out(1,:)=10000000000; %penalizing if dimensions don't match
else
out(1,:)=norm((vel_phase_fund - vel_phase_fund_for),2);
end

[status,cmdout]=unix ("forward_script2.scr"); %for first mode
vel_phase_first_for=load('forward_first_phase_velocity.txt');

if length(vel_phase_first_for)~= length(vel_phase_first)
out(2,:)=10000000000; %penalizing if dimensions don't match
else
out(2,:)=norm((vel_phase_first - vel_phase_first_for),2);
end

end
