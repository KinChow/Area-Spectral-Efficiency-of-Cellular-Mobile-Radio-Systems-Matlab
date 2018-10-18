clc;
clear;

% parameters
Ro=20/1000;  % Ro corresponds to the closest distance the mobile can be from the BS antenna, and is approximately equal to 20m for microcellular systems and 80m for macrocellular systems.
R_1=200/1000; % Cell Radius 200m
R_2=800/1000; % Cell Radius 800m
R_3=5000/1000; % Cell Radius 5000m 
R_4=100/1000:10^-2:1000/1000; % Cell Radius from 100m to 1000m
N_I_1=6;% The maximum number for n_i
N_I_2=2;
t=10^4; % repeat the process t time 
R_u_1=2:10^-1:10;
R_u_2=4;
R_u_3=6;
R_u_4=8;
a=2;
b_1=2;
b_2=4;
b_3=6;
h_B=10; % h_B=10 m for microcells and h_B=50 m for macrocells
h_m=2;
f_c_1=900*10^6; % 900MHz system and 2GHz system
lambda_c_1=3*10^8/f_c_1;
g_1=4*h_B*h_m/lambda_c_1/1000;
f_c_2=2000*10^6; % 900MHz system and 2GHz system
lambda_c_2=3*10^8/f_c_2;
g_2=4*h_B*h_m/lambda_c_2/1000;

% Comparison of the average uplink area spectral versus the normalized 
% reuse distance for different values of the additional path exponent b. 
[ASE_min_1,ASE_max_1,ASE_1]=uplink_1(Ro,R_1,N_I_1,t,R_u_1,a,b_1,g_1); % b=2
[ASE_min_2,ASE_max_2,ASE_2]=uplink_1(Ro,R_1,N_I_1,t,R_u_1,a,b_2,g_1); % b=4
[ASE_min_3,ASE_max_3,ASE_3]=uplink_1(Ro,R_1,N_I_1,t,R_u_1,a,b_3,g_1); % b=6
% plot
figure,
subplot(3,1,1), % b=2
plot(R_u_1,ASE_min_1,'b',R_u_1,ASE_max_1,'r',R_u_1,ASE_1,'g');
legend('Worst-Case Interference','Best-Case Interference','Simulations');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');
subplot(3,1,2), % b=4
plot(R_u_1,ASE_min_2,'b',R_u_1,ASE_max_2,'r',R_u_1,ASE_2,'g');
legend('Worst-Case Interference','Best-Case Interference','Simulations');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');
subplot(3,1,3), % b=6
plot(R_u_1,ASE_min_3,'b',R_u_1,ASE_max_3,'r',R_u_1,ASE_3,'g');
legend('Worst-Case Interference','Best-Case Interference','Simulations');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');

% Comparison of the average uplink area spectral versus the normalized 
% reuse distance for different cell sizes.
[ASE_min_4,ASE_max_4,ASE_4]=uplink_1(Ro,R_1,N_I_1,t,R_u_1,a,b_1,g_1); % R=200m
[ASE_min_5,ASE_max_5,ASE_5]=uplink_1(Ro,R_2,N_I_1,t,R_u_1,a,b_1,g_1); % R=800m
[ASE_min_6,ASE_max_6,ASE_6]=uplink_1(Ro,R_3,N_I_1,t,R_u_1,a,b_1,g_1); % b=5km
% plot
figure,
subplot(3,1,1), % R=200m
plot(R_u_1,ASE_min_4,'b',R_u_1,ASE_max_4,'r',R_u_1,ASE_4,'g');
legend('Worst-Case Interference','Best-Case Interference','Simulations');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');
subplot(3,1,2), % R=800m
plot(R_u_1,ASE_min_5,'b',R_u_1,ASE_max_5,'r',R_u_1,ASE_5,'g');
legend('Worst-Case Interference','Best-Case Interference','Simulations');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');
subplot(3,1,3), % R=5km
plot(R_u_1,ASE_min_6,'b',R_u_1,ASE_max_6,'r',R_u_1,ASE_6,'g');
legend('Worst-Case Interference','Best-Case Interference','Simulations');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');

% Average uplink area spectral effeciency versus cell radius for different
% reuse distances and carrier frequencies
[ASE_7]=uplink_2(Ro,R_4,N_I_1,t,R_u_2,a,b_1,g_1);
[ASE_8]=uplink_2(Ro,R_4,N_I_1,t,R_u_3,a,b_1,g_1);
[ASE_9]=uplink_2(Ro,R_4,N_I_1,t,R_u_4,a,b_1,g_1);
[ASE_10]=uplink_2(Ro,R_4,N_I_1,t,R_u_2,a,b_1,g_2);
[ASE_11]=uplink_2(Ro,R_4,N_I_1,t,R_u_3,a,b_1,g_2);
[ASE_12]=uplink_2(Ro,R_4,N_I_1,t,R_u_4,a,b_1,g_2);
figure,
semilogy(R_4,ASE_7,'b',R_4,ASE_10,'--b',R_4,ASE_8,'r',R_4,ASE_11,'--r',R_4,ASE_9,'g',R_4,ASE_12,'--g');
legend('fc=900 MHz','fc=2 GHz','fc=900 MHz','fc=2 GHz','fc=900 MHz','fc=2 GHz');
grid;
xlabel('Cell Radius R [Km]');
ylabel('Average Area Spectal Efficiency <ASE> [Bits/Sec/Hz/Km^2]');


% Average uplink area spectral efficiency versus normalized reuse distance
% with 120 cell sectorization
[ASE_min_13,ASE_max_13,ASE_13]=uplink_1(Ro,R_1,N_I_2,t,R_u_1,a,b_1,g_1); 
figure,
plot(R_u_1,ASE_min_13,'--b',R_u_1,ASE_min_1,'b',R_u_1,ASE_max_13,'--r',R_u_1,ASE_max_1,'r');
legend('3 Sectors','No Sectorization','3 Sectors','No Sectorization');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');

function [ASE_min,ASE_max,ASE]=uplink_1(Ro,R,N_I,t,R_u,a,b,g) 
% the position of desired user
u=rand(t,1);
r=Ro+(R-Ro)*sqrt(u);

% the polar coordinates (x,theta) of the N_I cochannel interferers 
x=zeros(t,N_I);
theta=zeros(t,N_I);
for n_i=1:N_I
    u_i=rand(t,1);
    v_i=rand(t,1);
    x(:,n_i)=Ro+(R-Ro)*sqrt(u_i);
    theta(:,n_i)=2*pi*v_i;
end

% the distance r_i from each cochannel interferer to the considered BS
D=R_u*R; % R_u normalized reuse distance
r_i=zeros(t,max(size(D)),N_I);
r_min=zeros(t,max(size(D)));
r_max=zeros(t,max(size(D)));
for n_i=1:N_I
    r_i(:,:,n_i)=sqrt(D.^2+x(:,n_i).^2+2.*D.*x(:,n_i).*sin(theta(:,n_i)));
end
for i=1:t
    for j=1:max(size(D))
        r_min(i,j)=D(j)-R;
        r_max(i,j)=D(j)+R;
    end
end

% the CIR of the desired user gamma_d
% gamma_d range from 1/N_I*(R*(R-1)/r)^a*((g+(R_u-1)*R)/(g+r))^b to 1/N_I*(R*(R+1)/r)^a*((g+(R_u+1)*R)/(g+r))^b
S_d=1./(r.^a.*(1+r/g).^b);
S_i=zeros(t,max(size(D)),N_I);
S_i_max=zeros(t,max(size(D)),N_I);
S_i_min=zeros(t,max(size(D)),N_I);
S_I=zeros(t,max(size(D)));
S_I_max=zeros(t,max(size(D)));
S_I_min=zeros(t,max(size(D)));
for n_i=1:N_I
    S_i(:,:,n_i)=1./(r_i(:,:,n_i).^a.*(1+r_i(:,:,n_i)/g).^b);
    S_i_max(:,:,n_i)=1./(r_max.^a.*(1+r_max/g).^b);
    S_i_min(:,:,n_i)=1./(r_min.^a.*(1+r_min/g).^b);
    S_I=S_I+S_i(:,:,n_i);
    S_I_max=S_I_max+S_i_max(:,:,n_i);
    S_I_min=S_I_min+S_i_min(:,:,n_i);
end
gamma_d=S_d./S_I;
gamma_d_max=S_d./S_I_max;
gamma_d_min=S_d./S_I_min;

% calculate the ASE
% ASE=mean(4/(pi*R_u*R)*log2(1+gamma_d));
ASE_min=mean(4./(pi.*D.^2).*log2(1+gamma_d_min));
ASE_max=mean(4./(pi.*D.^2).*log2(1+gamma_d_max));
ASE=mean(4./(pi.*D.^2).*log2(1+gamma_d));

end

function [ASE]=uplink_2(Ro,R,N_I,t,R_u,a,b,g) 
% the position of desired user
u=rand(t,1);
r=Ro+(R-Ro).*sqrt(u);

% the polar coordinates (x,theta) of the N_I cochannel interferers 
x=zeros(t,N_I);
theta=zeros(t,N_I);
for n_i=1:N_I
    for i=1:max(size(R))
        u_i=rand(t,1);
        v_i=rand(t,1);
        x(:,i,n_i)=Ro+(R(:,i)-Ro).*sqrt(u_i);
        theta(:,n_i)=2*pi*v_i;
    end
end

% the distance r_i from each cochannel interferer to the considered BS
D=R_u*R;
r_i=zeros(t,max(size(D)),N_I);
for n_i=1:N_I
    %for i=1:max(size(R))
        r_i(:,:,n_i)=sqrt(D.^2+x(:,:,n_i).^2+2.*D.*x(:,:,n_i).*sin(theta(:,n_i)));
    %end
end

% the CIR of the desired user gamma_d
% gamma_d range from 1/N_I*(R*(R-1)/r)^a*((g+(R_u-1)*R)/(g+r))^b to 1/N_I*(R*(R+1)/r)^a*((g+(R_u+1)*R)/(g+r))^b
S_d=1./(r.^a.*(1+r/g).^b);
S_I=zeros(t,max(size(D)));
for n_i=1:N_I
    S_I=S_I+1./(r_i(:,:,n_i).^a.*(1+r_i(:,:,n_i)/g).^b);
end
gamma_d=S_d./S_I;

% calculate the ASE
ASE=mean(4./(pi.*D.^2).*log2(1+gamma_d));
end