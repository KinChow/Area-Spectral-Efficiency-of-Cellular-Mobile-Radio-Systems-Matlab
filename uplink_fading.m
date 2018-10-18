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
R_u=2:0.1:10;
a=2;
b=2;
m_d=3; % m_d=1;m_d=2;m_d=3;
m_i=2; % m_i=1;m_i=2;m_i=3;
m_I=m_i;
t=10^5;
rho=m_d+m_I*N_I;

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
    for j=1:max(size(R_u))
        r_min(i,j)=D(j)-R;
        r_max(i,j)=D(j)+R;
    end
end

% the CIR of the desired user gamma_d
% gamma_d range from 1/N_I*(R*(R-1)/r)^a*((g+(R_u-1)*R)/(g+r))^b to 1/N_I*(R*(R+1)/r)^a*((g+(R_u+1)*R)/(g+r))^b
gamma_d_min=1/N_I.*(r_min./r).^a.*((g+r_min)./(g+r)).^b;
gamma_d_max=1/N_I.*(r_max./r).^a.*((g+r_max)./(g+r)).^b;
S_d=1./(r.^a.*(1+r/g).^b);
S_d_fading=gamrnd(m_d,S_d/m_d);

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
    S_i_fading(:,:,n_i)=gamrnd(m_i,S_i(:,:,n_i)/m_i);
    S_i_max_fading(:,:,n_i)=gamrnd(m_i,S_i_max(:,:,n_i)/m_i);
    S_i_min_fading(:,:,n_i)=gamrnd(m_i,S_i_min(:,:,n_i)/m_i);
    S_I_fading=S_I_fading+S_i_fading(:,:,n_i);
    S_I_max_fading=S_I_max_fading+S_i_max_fading(:,:,n_i);
    S_I_min_fading=S_I_min_fading+S_i_min_fading(:,:,n_i);
end
gamma_d=S_d./S_I;
gamma_d_fading=S_d_fading./S_I_fading;
gamma_d_max_fading=S_d_fading./S_I_max_fading;
gamma_d_min_fading=S_d_fading./S_I_min_fading;

% calculate the ASE
ASE_min=mean(4./(pi.*D.^2).*log2(1+gamma_d_min));
ASE_min_fading=mean(4./(pi.*D.^2).*log2(1+gamma_d_min_fading));
ASE_max=mean(4./(pi.*D.^2).*log2(1+gamma_d_max));
ASE_max_fading=mean(4./(pi.*D.^2).*log2(1+gamma_d_max_fading));
ASE=mean(4./(pi.*D.^2).*log2(1+gamma_d));
ASE_fading=mean(4./(pi.*D.^2).*log2(1+gamma_d_fading));

figure,
plot(R_u,ASE_min,'b',R_u,ASE_max,'r',R_u,ASE,'g',R_u,ASE_min_fading,'--b',R_u,ASE_max_fading,'--r',R_u,ASE_fading,'--g');
legend('Best-Case Interference No Fading','Simulations No Fading','Worst-Case Interference No Fading','Best-Case Interference Fading','Simulations Fading','Worst-Case Interference Fading')
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');