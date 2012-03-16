%% Protocol


N_scan_total = size(data{1},2);

for ii=1:N_scan_total
the_speeds(1,1,ii)=vit_alpha_mean(ii);
the_speeds(2,1,ii)=vit_beta_mean(ii);
end

unique_speed=zeros(2,4);

speeds_sign=sign(the_speeds);
unique_speed_sign=[1 -1 1 -1; 1 -1 -1 1];
unique_speed(1,1:2)=mean(abs(the_speeds(1,1,1:N_scan_total/2)));
unique_speed(2,1:2)=mean(abs(the_speeds(2,1,1:N_scan_total/2)));
unique_speed(1,3:4)=mean(abs(the_speeds(1,1,N_scan_total/2+1:N_scan_total)));
unique_speed(2,3:4)=mean(abs(the_speeds(2,1,N_scan_total/2+1:N_scan_total)));

unique_speed = unique_speed.*unique_speed_sign;
for ii=1:N_scan_total
    temp=speeds_sign(:,:,ii)'*unique_speed_sign;
    ind=find(temp==2);
    the_speeds(:,1,ii) = unique_speed(:,ind);
end

Nspeed = length(unique_speed(1,:))

