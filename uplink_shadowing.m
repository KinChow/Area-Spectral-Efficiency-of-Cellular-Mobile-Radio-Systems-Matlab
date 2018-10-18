clc;
clear;

% parameters
xi=10/log(10);
N_I=6;
R=200/1000;
Ro=20/1000;
h_B=10; % h_B=10 m for microcells and h_B=50 m for macrocells
h_m=2;
f_c=900*10^6; % 900MHz system and 2GHz system
lambda_c=3*10^8/f_c;
g=4*h_B*h_m/lambda_c/1000;
R_u_1=2:10^-1:10;
R_u_2=3:10^-1:10;
a=2;
b=2;
sigma_d=4;
sigma_I=4;
t=3*10^5;

D_1=R_u_1*R; % R_u for best
D_2=R_u_2*R; % R_u for worst
r_max=D_1+R;
r_min=D_2-R;
[ASE_max_shadowing_up,ASE_max_shadowing_low]=ulb(Ro,t,R,D_1,r_max,xi,sigma_d,sigma_I,N_I,a,b,g);
[ASE_min_shadowing_up,ASE_min_shadowing_low]=ulb(Ro,t,R,D_2,r_min,xi,sigma_d,sigma_I,N_I,a,b,g);

[ASE_min,ASE_max,ASE]=uplink_1(Ro,R,N_I,t,R_u_1,a,b,g);
[ASE_min_shadowing,ASE_max_shadowing,ASE_shadowing]=up_shadowing(Ro,R,N_I,xi,t,R_u_1,a,b,g,sigma_d,sigma_I);

%plot
figure,
plot(R_u_1,ASE_max_shadowing_up,'--r',R_u_2,ASE_min_shadowing_up,'--b');
legend('Best-Case Interference Upper Bound','Worst-Case Interference Upper Bound');
hold on,
plot(R_u_1,ASE_max_shadowing_low,'-.r',R_u_2,ASE_min_shadowing_low,'-.b');
legend('Best-Case Interference Lower Bound','Worst-Case Interference Lower Bound');
hold on,
plot(R_u_1,ASE_max_shadowing,'r',R_u_1,ASE_min_shadowing,'b',R_u_1,ASE_shadowing,'g');
legend('Best-Case Interference','Worst-Case Interference','Simulations');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');

figure,
plot(R_u_1,ASE_max,'r',R_u_1,ASE,'g',R_u_1,ASE_min,'b');
legend('Best-Case Interference No Shadowing','Simulations No Shadowing','Worst-Case Interference No Shadowing');
hold on,
plot(R_u_1,ASE_max_shadowing,'--r',R_u_1,ASE_shadowing,'--g',R_u_1,ASE_min_shadowing,'--b');
legend('Best-Case Interference Shadowing','Simulations Shadowing','Worst-Case Interference Shadowing');
grid;
xlabel('Normalized Reuse Distance Ru');
ylabel('<ASE> [Bits/Sec/Hz/Km^2]');

% ASE upper bound and lower bound
function [up,low]=ulb(Ro,t,R,D,r_m,xi,sigma_d,sigma_I,N_I,a,b,g)
u=rand(t,1);
r=Ro+(R-Ro)*sqrt(u);

sigma_S_I=sqrt(xi^2*log((N_I-1+exp(sigma_I^2/xi^2))/N_I));
sigma_gamma_d=sqrt(sigma_d^2+sigma_S_I^2);


mu_gamma_d_m=xi*log((r_m./r).^a.*((g+r_m)./(g+r)).^b)-xi*log(N_I)+(sigma_S_I^2-sigma_I^2)/(2*xi);

m=size(mu_gamma_d_m);

C_m_up=log2(exp(1))*(mu_gamma_d_m/xi+exp(-mu_gamma_d_m/xi+sigma_gamma_d^2/(2*xi^2)));

C_m_low=zeros(m(1),m(2));
for i=1:m(2)
    for j=1:m(1)
        C_m_low(j,i)=log2(exp(1))*(mu_gamma_d_m(j,i)./xi+Q(mu_gamma_d_m(j,i)./sigma_gamma_d)-exp(mu_gamma_d_m(j,i)/xi+sigma_gamma_d^2/(2*xi^2))*Q(mu_gamma_d_m(j,i)/sigma_gamma_d+sigma_gamma_d/xi));
    end
end

up=mean(4*C_m_up./(pi.*D.^2));
low=mean(4*C_m_low./(pi.*D.^2));
end

% ASE without shadowing
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

% ASE with shadowing
function [ASE_min_shadowing,ASE_max_shadowing,ASE_shadowing]=up_shadowing(Ro,R,N_I,xi,t,R_u,a,b,g,sigma_d,sigma_I)
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
S_d=1./(r.^a.*(1+r/g).^b);
mu_d=xi.*log(S_d);
S_d_shadowing=lognrnd(mu_d./xi,sigma_d./xi);
mu_I=zeros(t,max(size(D)),N_I);
mu_I_max=zeros(t,max(size(D)),N_I);
mu_I_min=zeros(t,max(size(D)),N_I);
S_i=zeros(t,max(size(D)),N_I);
S_i_max=zeros(t,max(size(D)),N_I);
S_i_min=zeros(t,max(size(D)),N_I);
S_I_shadowing=zeros(t,max(size(R_u)));
S_I_max_shadowing=zeros(t,max(size(D)));
S_I_min_shadowing=zeros(t,max(size(D)));
for n_i=1:N_I
    S_i(:,:,n_i)=1./(r_i(:,:,n_i).^a.*(1+r_i(:,:,n_i)/g).^b);
    mu_I(:,:,n_i)=xi.*log(S_i(:,:,n_i));
    S_i_max(:,:,n_i)=1./(r_max.^a.*(1+r_max/g).^b);
    mu_I_max(:,:,n_i)=xi.*log(S_i_max(:,:,n_i));
    S_i_min(:,:,n_i)=1./(r_min.^a.*(1+r_min/g).^b);
    mu_I_min(:,:,n_i)=xi.*log(S_i_min(:,:,n_i));
    S_I_shadowing=S_I_shadowing+lognrnd(mu_I(:,:,n_i)./xi,sigma_I./xi);
    S_I_max_shadowing=S_I_max_shadowing+lognrnd(mu_I_max(:,:,n_i)./xi,sigma_I./xi);
    S_I_min_shadowing=S_I_min_shadowing+lognrnd(mu_I_min(:,:,n_i)./xi,sigma_I./xi);
end

gamma_d_shadowing=S_d_shadowing./S_I_shadowing;
gamma_d_max_shadowing=S_d_shadowing./S_I_max_shadowing;
gamma_d_min_shadowing=S_d_shadowing./S_I_min_shadowing;
    
% calculate the ASE
ASE_shadowing=mean(4./(pi.*D.^2).*log2(1+gamma_d_shadowing));
ASE_max_shadowing=mean(4./(pi.*D.^2).*log2(1+gamma_d_max_shadowing));
ASE_min_shadowing=mean(4./(pi.*D.^2).*log2(1+gamma_d_min_shadowing));
end

% Q function
function Qz=Q(z)
fun = @(x) exp(-x.^2./2);
Qz=1/sqrt(2*pi)*integral(fun,z,Inf);
end

% Desired User Average Capacity
function C=DUAC(mu_gamma_d,sigma_gamma_d,xi)
A=log2(exp(1))*xi/(sqrt(2*pi)*sigma_gamma_d);
fun = @(x) log(1+x)./x.*exp(-(xi.*log(x)-mu_gamma_d).^2/(2*sigma_gamma_d^2));
B=integral(fun,0,Inf);
C=A*B;
end