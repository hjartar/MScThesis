clear all;
close all;
tic
cd /Users/hjorturhjartar/Documents/POLIT/Thesis/Matlab
data_dir = '/Users/hjorturhjartar/Documents/POLIT/Thesis/Matlab'

%variables
global eta rho w0 a r A rl  rd liq

% deposits
rho = 0.5;
a = 0.5; 
w0 = 1;
liq = 1

%loans
A=0.5;



%options
options = optimset('TolX',0.0000000000000000000000000000000001,'Display','off');
format long

%===========
%Deposits market
%===========

%matrix
etaMat= [0.5 1.5 5 20];
[K,k] = size(etaMat);
rdDiff = 0.01;
rdMat= -0.5:rdDiff:0.5;
rdMat = rdMat';
[N,n] = size(rdMat);
detaMat = zeros(N,k);
dMat = zeros(N,1);

%deposits for elasticity calc
for j = 1:k
    eta = etaMat(1,j);
   
    for i =1:N
        rd = rdMat(i,1);
        dMat(i,1) = fminbnd(@deposits1,0,w0,options);
        detaMat(i,j) = dMat(i,1);       
    end    
end

%deposit elasticity
logdMat = log(detaMat);
ElasticityMat = diff(logdMat)./rdDiff;

%Implied risk free rate
rfMat = rdMat(2:N,1) - 1./ElasticityMat*0.01;

% plot(rdMat,logdMat);
plot(rdMat, detaMat);
grid on


%===========
%Wealth dependent deposits
%===========%===========



%variables
global eta rho w0 a r A rl  rd liq

% deposits
rho = 0.5;
a = 0.5; 

liq = 1

%loans
A=0.5;

%options
options = optimset('TolX',0.0000000000000000000000000000000001,'Display','off');
format long

%matrix
etaMat= [20]%0.1:0.2:2;
wMat= [1 1.2 1.5 2]%0.1:0.2:2;
[K,k] = size(wMat);
rdDiff = 0.01;
rdMat= -0.5:rdDiff:0.5;
rdMat = rdMat';
[N,n] = size(rdMat);
detaMat1 = zeros(N,k);
dMat = zeros(N,1);

%deposits for elasticity calc
for j = 1:k
    w0 = wMat(1,j);
   eta=etaMat;
    for i =1:N
        rd = rdMat(i,1);
        dMat(i,1) = fminbnd(@deposits2,0,w0,options);
        detaMat1(i,j) = dMat(i,1)/w0;       
    end    
end

%deposit elasticity
logdMat = log(detaMat1);
ElasticityMat = diff(logdMat)./rdDiff;

%Implied risk free rate
rfMat = rdMat(2:N,1) - 1./ElasticityMat*0.01;

% plot(rdMat,logdMat);
plot(rdMat, detaMat1);
grid on

set(gca,'FontSize',30)


f3 = figure(1);


subplot(1,2,1)
plot(rdMat, detaMat);
legend('\eta = 0.5','\eta = 1.5','\eta = 5','\eta = 20','Location','southeast') %[0.5 1.5 5 20]
grid on
xlabel('Net deposit rate: R_D-1')
ylabel('Share of wealth in deposits')
title('Deposit allocation at varying elasticities of substitution')
set(gca,'FontSize',16)
set(gca,'xtick',[-0.5:0.1:0.5])
pbaspect([1 1 1])

subplot(1,2,2)
plot(rdMat, detaMat1);
legend('w_0 = 1','w_0 = 1.2','w_0 = 1.5','w_0 = 2','Location','southeast') %[1 1.2 1.5 2]
grid on
xlabel('Net deposit rate: R_D-1')
ylabel('Share of wealth in deposits')
title('Deposit allocation at varying amounts of initial wealth')
set(gca,'FontSize',16)
set(gca,'xtick',[-0.5:0.1:0.5])
pbaspect([1 1 1])
set(gca,'FontSize',16)


saveas(f3, strcat(data_dir, '\deposits.eps'),'epsc')




