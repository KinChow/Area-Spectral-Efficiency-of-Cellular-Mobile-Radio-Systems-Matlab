clc;
clear;

% parameters
N_I=6;
R=200/1000;
Ro=20/1000;
h_B=10; % h_B=10 m for microcells and h_B=50 m for macrocells
h_m=2;
f_c=900*10^6; % 900MHz system and 2GHz system
lambda_c=3*10^8/f_c;
g=4*h_B*h_m/lambda_c/1000;
R_u=2:0.5:10;
a=2;
b=2;
m_d=3;
m_i=1;
m_I=m_i;
t=3*10^5;
rho=m_d+m_I*N_I;
sigma_d=6;
sigma_I=6;
xi=10/log(10);

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
%gamma_d_min=1/N_I.*(r_min./r).^a.*((g+r_min)./(g+r)).^b;
%gamma_d_max=1/N_I.*(r_max./r).^a.*((g+r_max)./(g+r)).^b;
S_d=1./(r.^a.*(1+r/g).^b);

mu_d=xi*log(S_d);
mu_d_combined=xi*(psi(m_d)-log(m_d))+mu_d;

sigma_d_combined=sqrt(xi^2*rz(2,m_d)+sigma_d^2);

S_d_combined=lognrnd(mu_d_combined./xi,sigma_d_combined./xi);

S_I=zeros(t,max(size(D)));
S_I_combined=zeros(t,max(size(D)));
S_I_max_combined=zeros(t,max(size(D)));
S_I_min_combined=zeros(t,max(size(D)));
for n_i=1:N_I
    S_i(:,:,n_i)=1./(r_i(:,:,n_i).^a.*(1+r_i(:,:,n_i)/g).^b);
    S_i_max=1./(r_max.^a.*(1+r_max/g).^b);
    S_i_min=1./(r_min.^a.*(1+r_min/g).^b);
    
    mu_I(:,:,n_i)=xi.*log(S_i(:,:,n_i));
    mu_I_max=xi.*log(S_i_max);
    mu_I_min=xi.*log(S_i_min);
    
    mu_I_combined(:,:,n_i)=xi*(psi(m_I)-log(m_I))+mu_I(:,:,n_i);
    mu_I_max_combined=xi*(psi(m_I)-log(m_I))+mu_I_max;
    mu_I_min_combined=xi*(psi(m_I)-log(m_I))+mu_I_min;
    
    sigma_I_combined=sqrt(xi^2*rz(2,m_I)+sigma_I^2);
    
    S_i_combined(:,:,n_i)=lognrnd(mu_I_combined(:,:,n_i)./xi,sigma_I_combined./xi);
    S_i_max_combined(:,:,n_i)=lognrnd(mu_I_max_combined./xi,sigma_I_combined./xi);
    S_i_min_combined(:,:,n_i)=lognrnd(mu_I_min_combined./xi,sigma_I_combined./xi);
    
    S_I=S_I+S_i(:,:,n_i);
    S_I_combined=S_I_combined+S_i_combined(:,:,n_i);
    S_I_max_combined=S_I_max_combined+S_i_max_combined(:,:,n_i);
    S_I_min_combined=S_I_min_combined+S_i_min_combined(:,:,n_i);
end

gamma_d=S_d./S_I;
gamma_d_combined=S_d_combined./S_I_combined;
gamma_d_max_combined=S_d_combined./S_I_max_combined;
gamma_d_min_combined=S_d_combined./S_I_min_combined;

[ASE_min,ASE_max,ASE]=montecarlo1(Ro,R,N_I,t,R_u,a,b,g); 
ASE_combined=mean(4./(pi.*D.^2).*log2(1+gamma_d_combined));
ASE_max_combined=mean(4./(pi.*D.^2).*log2(1+gamma_d_max_combined));
ASE_min_combined=mean(4./(pi.*D.^2).*log2(1+gamma_d_min_combined));

figure,
plot(R_u,ASE,'g',R_u,ASE_combined,'--g');
hold on,
plot(R_u,ASE_max,'r',R_u,ASE_max_combined,'--r');
hold on,
plot(R_u,ASE_min,'b',R_u,ASE_min_combined,'--b');
legend('No Shadowing/Fading','Shadowing and Fading');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');

function result=rz(z,n)
if n==1
    result=zeta(z);
elseif n>1
    s=0;
    for i=1:n
        s=s+1/i^z;
    end
    result=zeta(z)-s;
else
    warning('cannot compute value')
end
end

function [ASE_min,ASE_max,ASE]=montecarlo1(Ro,R,N_I,t,R_u,a,b,g) 
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
r_min=D-R;
r_max=D+R;
for n_i=1:N_I
    r_i(:,:,n_i)=sqrt(D.^2+x(:,n_i).^2+2.*D.*x(:,n_i).*sin(theta(:,n_i)));
end

% the CIR of the desired user gamma_d
% gamma_d range from 1/N_I*(R*(R-1)/r)^a*((g+(R_u-1)*R)/(g+r))^b to 1/N_I*(R*(R+1)/r)^a*((g+(R_u+1)*R)/(g+r))^b
gamma_d_min=1/N_I.*(r_min./r).^a.*((g+r_min)./(g+r)).^b;
gamma_d_max=1/N_I.*(r_max./r).^a.*((g+r_max)./(g+r)).^b;
S_d=1./(r.^a.*(1+r/g).^b);
S_I=zeros(t,max(size(R_u)));
for n_i=1:N_I
    S_I=S_I+1./(r_i(:,:,n_i).^a.*(1+r_i(:,:,n_i)/g).^b);
end
gamma_d=S_d./S_I;

% calculate the ASE
% ASE=mean(4/(pi*R_u*R)*log2(1+gamma_d));
ASE_min=mean(4./(pi.*D.^2).*log2(1+gamma_d_min));
ASE_max=mean(4./(pi.*D.^2).*log2(1+gamma_d_max));
ASE=mean(4./(pi.*D.^2).*log2(1+gamma_d));

end