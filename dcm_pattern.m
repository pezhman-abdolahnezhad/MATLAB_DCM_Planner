clc
clear
rf=[0 0.115 0; 0 -0.115 0;.5 .115 0;1 -.115 0;1.5 .115 0;2 -.115 0;2.5 .115 0];
delta_z_vrp = 0.8;
g = 9.81;
t_step = .8;
t_ds = .3;
alpha = .3;
% for i=1:6
figure;
plot(rf(:,1),rf(:,2),'.');
xlabel( 'X');
ylabel('Y');
grid on;
r_vrpd=rf;
r_vrpd =[r_vrpd(:,1),r_vrpd(:,2), r_vrpd(:,3) + delta_z_vrp];
figure
scatter3(r_vrpd(:,1),r_vrpd(:,2), r_vrpd(:,3),'r','filled');
xlabel( 'X');
ylabel('Y');
zlabel('Z');
grid on;
hold on;
scatter3(rf(:,1),rf(:,2),rf(:,3),'g','filled');
hold on
%% desired DCM location
ksi_d_eos = zeros(7,3);
ksi_d_eos(7,:) =  r_vrpd (7,:);
for i=6:-1:1
    ksi_d_eos(i,:) = r_vrpd(i+1,:) + (exp((-sqrt(g/delta_z_vrp)*t_step)))*(ksi_d_eos(i+1,:)-r_vrpd(i+1,:));
end 
scatter3(ksi_d_eos(:,1),ksi_d_eos(:,2),ksi_d_eos(:,3),'b')
hold on
% ;kplot3(r_vrpd(:,1),r_vrpd(:,2), r_vrpd(:,3),'r')
sample_time = 0.005;
t = 0:sample_time:t_step * 7 - sample_time;
ksi_d = zeros(length(t), 3);

for j = 1:length(t)
    step_num = ceil(j*sample_time/t_step);
    local_t=mod(t(j),t_step);
    ksi_d(j,:) = r_vrpd(step_num,:) + exp(sqrt(g / delta_z_vrp)*(local_t - t_step)).*(ksi_d_eos(step_num, :) - r_vrpd(step_num, :));
    
%   hold on
end
plot3(ksi_d(:,1),ksi_d(:,2),ksi_d(:,3))
% scatter3(ksi_d(:,1),ksi_d(:,2),ksi_d(:,3),'b','filled');
%% double support
% ksi_iniDS = r_vrpd + exp(-sqrt(g/delta_z_vrp)*t_ds)*(ksi_d-r_vrpd)
ksi_iniDS(1,:) = ksi_d(1,:)
ksi_eoDS(i,:) = r_vrpd(i,:) + exp(sqrt(g/delta_z_vrp)*t_ds*(1-alpha))*(ksi_d(i,:)-r_vrpd(i,:))
for i = 1:7
    if i==1
        ksi_iniDS(i,:) = ksi_d(i,:)
        ksi_eoDS(i,:) = r_vrpd(i,:) + exp(sqrt(g/delta_z_vrp)*t_ds*(1-alpha))*(ksi_d(i,:)-r_vrpd(i,:))
    else
        ksi_eoDS(i,:) = r_vrpd(i,:) + exp(sqrt(g/delta_z_vrp)*t_ds*(1-alpha))*(ksi_d_eos(i-1,:)-r_vrpd(i,:))
        ksi_iniDS(i,:) = r_vrpd(i-1,:) + exp(-sqrt(g/delta_z_vrp)*t_ds*alpha)*(ksi_d_eos(i-1,:)-r_vrpd(i-1,:));
    end
end
scatter3(ksi_eoDS(:,1),ksi_eoDS(:,2),ksi_eoDS(:,3),'y','filled');
hold on
scatter3(ksi_iniDS(:,1),ksi_iniDS(:,2),ksi_iniDS(:,3),'k','filled');
%% polynomial 
for i= 1:7
    if i==1
        ksi_dot_iniDS(i,:) = (ksi_iniDS(i,:) - r_vrpd(i,:))*sqrt(g/delta_z_vrp);
        ksi_dot_eoDS(i,:) = (ksi_eoDS(i,:) - r_vrpd(i,:))*sqrt(g/delta_z_vrp);
    else
        ksi_dot_iniDS(i,:) = (ksi_iniDS(i,:) - r_vrpd(i-1,:))*sqrt(g/delta_z_vrp);
        ksi_dot_eoDS(i,:) = (ksi_eoDS(i,:) - r_vrpd(i,:))*sqrt(g/delta_z_vrp);
    end
end
%%

ksi_traj_ds = struct('index', {}, 'ksi_ds', {});
for i=1:7
    wpts=[ksi_iniDS(i,1) ksi_eoDS(i,1);ksi_iniDS(i,2) ksi_eoDS(i,2)];
    velpts = [ksi_dot_iniDS(i,1) ksi_dot_eoDS(i,1);ksi_dot_iniDS(i,2) ksi_dot_eoDS(i,2)];
    if i == 1

        tpts=[0 (1-alpha) * t_ds];
        tvec=tpts(1):sample_time:tpts(end);
        ksi_traj_ds(i).index = i;
        [q, qd, qdd, pp] = cubicpolytraj(wpts, tpts, tvec,'VelocityBoundaryCondition',velpts);
        ksi_traj_ds(i).ksi_ds = q';
    else
        
        tpts=[0 t_ds];
        tvec=tpts(1):sample_time:tpts(end);
        ksi_traj_ds(i).index = i;
        [q, qd, qdd, pp] = cubicpolytraj(wpts, tpts, tvec,'VelocityBoundaryCondition',velpts);
        ksi_traj_ds(i).ksi_ds = q';
    end
    plot(q(1, :), q(2, :))
end 
%% merge
final_ksi_traj = ksi_d;
final_ksi_traj(1:(t_ds * (1-alpha)) / sample_time + 1, 1:2) = ksi_traj_ds(1).ksi_ds;
for i=1:6
    initial_ds = floor((i * t_step - alpha * t_ds) / sample_time);
    final_ds = floor((i * t_step + (1 - alpha) * t_ds) / sample_time);
    final_ksi_traj(initial_ds:final_ds, 1:2) = ksi_traj_ds(i+1).ksi_ds;
end
% figure()
plot(final_ksi_traj(:, 1), final_ksi_traj(:, 2))
hold on

%% com pattern
t_s = 0:sample_time:t_step * 7 - sample_time;
com_traj_final(1,:)=[0,0,.8]
for index = 2:length(t)
    xx = zeros(1,3);
    for t =1:index
        inte =1/100*final_ksi_traj(t,:)* exp((t/100)/sqrt(delta_z_vrp/g));
        xx = inte + xx;
    end
    com_traj_final(index,:) = (xx/(sqrt(delta_z_vrp/g)) + com_traj_final(1,:))*exp((-index/100)/sqrt(delta_z_vrp/g));  
end
  plot(com_traj_final(:,1),com_traj_final(:,2));
    
