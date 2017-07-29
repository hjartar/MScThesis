clear all;
close all;
tic
format long g
cd /Users/hjorturhjartar/Documents/POLIT/Thesis/Matlab
data_dir = '/Users/hjorturhjartar/Documents/POLIT/Thesis/Matlab'

%variables
global eta rho w0 a A liq rd rf lambda2solve

% deposits
eta = 50
rho = 20;
a = 0.5; 
w0 =20;
ws = w0;
eps = 0.01;
liq = 1;
%loans
gamma = 0.5;
A=0.2;
eL= 2;
E0=0.2;
reps = 2;
lcost =0

int=[-0.1 0.1]
%options
options = optimset('TolX',0.0000000000000000000000000000000001,'Display','off');
format long g

rdDiff = 0.01;
rdMat= 0.5:rdDiff:1.5;
rdMat = rdMat';
[N,n] = size(rdMat);
detaMat = zeros(N,1);
dMat = zeros(N,4);
dElaMat = zeros(N,1);
dDevMat = zeros(N,1);
dElaAltMat = zeros(N,1);
cdMat = zeros(N,1);

%deposits for elasticity calc

    for i =1:N
        rd = rdMat(i,1);
        dMat(i,1) = fminbnd(@deposits1,0,w0,options);
        dMat(i,2) = log(dMat(i,1));
        rd= rdMat(i,1)-eps;
        dMat(i,3) = fminbnd(@deposits1,0,w0,options);
        dMat(i,4) = log(dMat(i,3));
        dDevMat(i,1) = (dMat(i,1)-dMat(i,3))/(eps);
        dElaMat(i,1) = (dMat(i,2)-dMat(i,4))/((eps)*dMat(i,1))*rdMat(i,1);
        dElaAltMat(i,1) = dDevMat(i,1)/dMat(i,1)*rdMat(i,1);
        detaMat(i,1) = dMat(i,1); 
    end    

%plot(detaMat)
%plot(mark)
mark = (dElaMat(:,1)./(dElaMat(:,1)+1))



rf_and_lambda = rdMat./mark;
[rdMat, rf_and_lambda]

dMarket = [rdMat, dMat(:,1), dElaMat(:,1), mark, rf_and_lambda];
dMarket = dMarket(~any(isnan(dMarket),2),:)

Drd = @(rd) interp1(dMarket(:,1),dMarket(:,2),rd,'spline');
eDrd = @(rd) interp1(dMarket(:,1),dMarket(:,3),rd,'spline');
cDrd = @(rd) interp1(dMarket(:,1),dMarket(:,2),rd,'spline')
rdD = @(d) interp1(dMarket(:,2),dMarket(:,1),d,'spline')

rd_rf = @(d) interp1(dMarket(:,5),dMarket(:,1),d,'spline')

rd_rf(1)


%fplot(cDrd, int)
%fplot(lambda1, [-0.5 0.1])

%plot(dMarket)
%rd_solve = @(rd) (eDrd(rd)/(eDrd(rd)+1))*(rf+lambda2solve)-rd;

%rf=0.05
%fzero(rd_solve,0.1,options)


%================
% Bank model






format short g
global lambda2solve
options = optimset('TolX',0.0000000000000000000000000000000001,'Display','off');
rfv =-0.05:0.005:0.15;
Ev = [5 2 1.5 1 0.5 0]

sol1 = zeros(length(rfv),12);
solL = zeros(length(rfv),length(Ev));
solM = zeros(length(rfv),length(Ev));
solE = zeros(length(rfv),length(Ev));
lcost=0.05
for k = 1:length(Ev)
E0=Ev(k)
for j = 1:length(rfv)
        rf=rfv(j);
    lambda2=0;
    
    lambda1 = @(E) ((E/gamma)^(A-1)*A/(eL/(eL-1))-rf-lambda2-lcost)/gamma;
    rl1 = @(E) max((eL/(eL-1))*(rf+max(lambda1(E),0)*gamma+lambda2+lcost),0);
    Lrl1 = @(E) (rl1(E)/A)^(1/(A-1));

              
    rd1 = rd_rf(rf+1);           
    Drd1 = Drd(rd1);
    
    solve1= @(E) abs(Lrl1(E)*(rl1(E)-rf-lcost)+Drd1*(1+rf-rd1) + E0*(1+rf) -E);

    Esolve = fminbnd(solve1,0,10,options);
    M = E0 + Drd1-Lrl1(Esolve);
    method = 0;
    Lnim = (rl1(Esolve)-rf) * Lrl1(Esolve);
    Dnim = Drd1 * ((rf+1)-rd1 );
   
        sol1(j,:) = [rf Esolve rl1(Esolve) Lrl1(Esolve) rd1 Drd1 M solve1(Esolve) 0 method Lnim Dnim];

   
end
solL(:,k)=sol1(:,4)
solM(:,k)=sol1(:,7)
solE(:,k)=sol1(:,2)-E0
end
plot(sol1(:,1),sol1(:,4),sol1(:,1),sol1(:,6))
%plot(sol1(:,1),sol1(:,2),sol1(:,1),sol1(:,7))
f3 = figure(1);
subplot(1,1,1)
plot(sol1(:,1),solL)
legend('E_0=5','E_0=2','E_0=1.5','E_0=1','E_0=0.5','E_0=0')
grid on
ylim([0 6])
pbaspect([1 1 1])
set(gca,'FontSize',14)
ylabel('Loans')
xlabel('R_M-1')
title('Lending under capital constraint')
saveas(f3, strcat(data_dir, '\bank12.eps'),'epsc')
f3 = figure(1);




