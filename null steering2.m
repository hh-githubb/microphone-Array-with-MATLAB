close all; clc; clear all;
hcfontsize = 20;
c=340;
theta_d=0/180*pi;
theta1=35/180*pi;
theta2=47/180*pi;
delta=1e-2;
alpha_values=[1e-3 1e-4 1e-5];
M_values=[4 6 10 12]; % Number of sensors
f=linspace(0.001,8e3,12)';  % frequency
W1_dB_values=zeros(length(f),length(M_values));
W2_dB_values=zeros(length(f),length(M_values));
W3_dB_values=zeros(length(f),length(M_values));
W4_dB_values=zeros(length(f),length(M_values));
i_c=[1 0 0]';
for idxM=1:length(M_values)
    M=M_values(idxM);
    [m,n]=meshgrid(1:M,1:M);
    for i=1:length(f)
        fk=f(i);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        C=[exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d)) exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta1)) exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta2))];
        Gamma0=sinc(2*fk*delta/c*(m-n));
        alpha=alpha_values(1);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;
        W1_dB_values(i,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        alpha=alpha_values(2);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;
        W2_dB_values(i,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        alpha=alpha_values(3);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;
        W3_dB_values(i,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        h=d/M;
        W4_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
    end
end
figure
plot(f/1e3,W1_dB_values(:,1),'-bo');
hold on
plot(f/1e3,W2_dB_values(:,1),'--g*');
plot(f/1e3,W3_dB_values(:,1),':rs');
plot(f/1e3,W4_dB_values(:,1),'-.c^');
hold off
xlabel('F');
ylabel('WNG');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(f/1e3,W1_dB_values(:,2),'-bo');
hold on
plot(f/1e3,W2_dB_values(:,2),'--g*');
plot(f/1e3,W3_dB_values(:,2),':rs');
plot(f/1e3,W4_dB_values(:,2),'-.c^');
hold off
xlabel('F');
ylabel('WNG');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(f/1e3,W1_dB_values(:,3),'-bo');
hold on
plot(f/1e3,W2_dB_values(:,3),'--g*');
plot(f/1e3,W3_dB_values(:,3),':rs');
plot(f/1e3,W4_dB_values(:,3),'-.c^');
hold off
xlabel('F');
ylabel('WNG');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(f/1e3,W1_dB_values(:,4),'-bo');
hold on
plot(f/1e3,W2_dB_values(:,4),'--g*');
plot(f/1e3,W3_dB_values(:,4),':rs');
plot(f/1e3,W4_dB_values(:,4),'-.c^');
hold off
xlabel('F');
ylabel('WNG');

