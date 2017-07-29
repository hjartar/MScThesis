clear all;
close all;
tic
format long g
cd /Users/hjorturhjartar/Documents/POLIT/Thesis/Matlab
data_dir = '/Users/hjorturhjartar/Documents/POLIT/Thesis/Matlab'

%variables
global eta rho w0 a A liq rd rf lambda2solve

% deposits
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

CoE =0.2


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



%================
% Bank model
%================

format long g
global lambda2solve
options = optimset('TolX',0.0000000000000000000000000000000001,'Display','off');

gamma
%interest rate environment
rf=0.05
t= 60
rfv = ones(1,t)*rf;
rs=10
re=50
rfv(rs:re)=+0.02

%rhi=0.04
%rlo=0.01
%rfv(rs:re)=rlo:(rhi-rlo)/(re-rs):rhi


sol = zeros(length(rfv),12);

%terminal bank values

rf_term =rf
E0_bank =0
rl_term = (eL/(eL-1))*(rf_term+lcost)
L_term =  (rl_term/A)^(1/(A-1))

    rd_term = rd_rf(rf+1);           
    D_term = Drd(rd_term);

E1_term = E0_bank*(1+rf_term) + L_term*(rl_term-rf) + D_term*(rf+1-rd_term)
div_pol_term = 1-E0_bank/E1_term
div_term = div_pol_term*E1_term
V_term = div_term/CoE+div_term

constraint1_check = (L_term*gamma<V_term)
constraint2_check = (- L_term + E0_bank +D_term>0)

% backward calculation

%vectors
V_vector = zeros(t,1)
V_vector(t,1) = V_term

div_vector = zeros(t,1)
div_vector(t,1) = V_term

E_vector = zeros(t,1)
E_vector(1,1) = E0_bank

M_vector = zeros(t,1)
L_vector = zeros(t,1)
Lstar_vector = zeros(t,1)

D_vector = zeros(t,1)
Lnim_vector = zeros(t,1)
Dnim_vector = zeros(t,1)
roc_vector = zeros(t,1)


B = zeros (t,1);

for i = 1:t
    B(i)= 1/(1+CoE)^(i-1);  %discount factor
end


for k = 1:2
for j = 1:t-1
        bV_t1 = B(2:t-j+1,1)'*div_vector(1+j:t,1);

    lambda2=0;
  rf=rfv(j);
    if j==1
        E0 = E0_bank;
         
    else
        E0 = E_vector(j);
    end
    lambda1 = @(V) ((V/gamma)^(A-1)*A/(eL/(eL-1))-rf-lambda2-lcost)/gamma;
    rl1 = @(V) (eL/(eL-1))*(rf+max(lambda1(V),0)*gamma+lambda2+lcost);
    Lrl1 = @(V) (rl1(V)/A)^(1/(A-1));
    Lstar_vector(j) =   ((eL/(eL-1))*(rf+lcost)/A)^(1/(A-1)) ;
    rd1 =  rd_rf(rf+1+lambda2);
    Drd1 =  Drd(rd1);
    
    Nim1 = @(V) Lrl1(V)*(rl1(V)-rf)+Drd1*(rf+1-rd1);
    
    E_Book = @(V) Nim1(V) + E_vector(j)*(1+rf);
    div_pol = @(V) max(1-E0_bank/E_Book(V),0);
    
    solve1= @(V)  abs(E_Book(V)*div_pol(V) + bV_t1 -V);

    %solve1= @(V)  E_Book(V)*div_pol(V) + bV_t1 -V;
    Vsolve = fminbnd(solve1,0,V_term,options);
    %Vsolve = fzero(solve1,E0_bank,options);
    div_solve =  E_Book(Vsolve)*div_pol(Vsolve);
    
    M = E0 + Drd1-Lrl1(Vsolve);
    method = 0;
    
    
   
        div_vector(j)=div_solve;
        V_vector(j)=Vsolve;
        M_vector(j) = M;
        L_vector(j) = Lrl1(Vsolve);
        D_vector(j) = Drd1;
        Lnim_vector(j) = Lrl1(Vsolve)*(rl1(Vsolve)-rf);
        Dnim_vector(j) = Drd1*(rf-rd1);

        E_vector(j+1)= E_Book(Vsolve) - div_solve;

end
end

toc
f3 = figure(1);

subplot(1,2,2)
hold on 
plot(1:t, V_vector)

plot(1:t, L_vector)
plot(1:t, Lstar_vector,'--')
pbaspect([1 1 1])

legend('Present value of bank','Constrained optimal lending','Unconstrained optimal lending','Location','southeast')


x = [0.75 0.75];
y = [0.6 0.52];
p= annotation('textarrow',x,y,'String','Lending increase is delayed ')
p.FontSize = 7;
xlabel('Time Period')
set(gca,'FontSize',7)
title('Long low interest rate period')
grid on
%a.FontSize = 8;
hold off



%=======================



rf=0.05
t= 60
rfv = ones(1,t)*rf;
rs=10
re=30
rfv(rs:re)=+0.02

%rhi=0.04
%rlo=0.01
%rfv(rs:re)=rlo:(rhi-rlo)/(re-rs):rhi


sol = zeros(length(rfv),12);

%terminal bank values

rf_term =rf
E0_bank =0
rl_term = (eL/(eL-1))*(rf_term+lcost)
L_term =  (rl_term/A)^(1/(A-1))

    rd_term = rd_rf(rf+1);           
    D_term = Drd(rd_term);

E1_term = E0_bank*(1+rf_term) + L_term*(rl_term-rf) + D_term*(rf+1-rd_term)
div_pol_term = 1-E0_bank/E1_term
div_term = div_pol_term*E1_term
V_term = div_term/CoE+div_term

constraint1_check = (L_term*gamma<V_term)
constraint2_check = (- L_term + E0_bank +D_term>0)

% backward calculation

%vectors
V_vector = zeros(t,1)
V_vector(t,1) = V_term

div_vector = zeros(t,1)
div_vector(t,1) = V_term

E_vector = zeros(t,1)
E_vector(1,1) = E0_bank

M_vector = zeros(t,1)
L_vector = zeros(t,1)
Lstar_vector = zeros(t,1)

D_vector = zeros(t,1)
Lnim_vector = zeros(t,1)
Dnim_vector = zeros(t,1)
roc_vector = zeros(t,1)


B = zeros (t,1);

for i = 1:t
    B(i)= 1/(1+CoE)^(i-1);  %discount factor
end


for k = 1:2
for j = 1:t-1
        bV_t1 = B(2:t-j+1,1)'*div_vector(1+j:t,1);

    lambda2=0;
  rf=rfv(j);
    if j==1
        E0 = E0_bank;
         
    else
        E0 = E_vector(j);
    end
    lambda1 = @(V) ((V/gamma)^(A-1)*A/(eL/(eL-1))-rf-lambda2-lcost)/gamma;
    rl1 = @(V) (eL/(eL-1))*(rf+max(lambda1(V),0)*gamma+lambda2+lcost);
    Lrl1 = @(V) (rl1(V)/A)^(1/(A-1));
    Lstar_vector(j) =   ((eL/(eL-1))*(rf+lcost)/A)^(1/(A-1)) ;
    rd1 =  rd_rf(rf+1+lambda2);
    Drd1 =  Drd(rd1);
    
    Nim1 = @(V) Lrl1(V)*(rl1(V)-rf)+Drd1*(rf+1-rd1);
    
    E_Book = @(V) Nim1(V) + E_vector(j)*(1+rf);
    div_pol = @(V) max(1-E0_bank/E_Book(V),0);
    
    solve1= @(V)  abs(E_Book(V)*div_pol(V) + bV_t1 -V);

    %solve1= @(V)  E_Book(V)*div_pol(V) + bV_t1 -V;
    Vsolve = fminbnd(solve1,0,V_term,options);
    %Vsolve = fzero(solve1,E0_bank,options);
    div_solve =  E_Book(Vsolve)*div_pol(Vsolve);
    
    M = E0 + Drd1-Lrl1(Vsolve);
    method = 0;
    
    
   
        div_vector(j)=div_solve;
        V_vector(j)=Vsolve;
        M_vector(j) = M;
        L_vector(j) = Lrl1(Vsolve);
        D_vector(j) = Drd1;
        Lnim_vector(j) = Lrl1(Vsolve)*(rl1(Vsolve)-rf);
        Dnim_vector(j) = Drd1*(rf-rd1);

        E_vector(j+1)= E_Book(Vsolve) - div_solve;

end
end

toc


subplot(1,2,1)
hold on 
plot(1:t, V_vector)
plot(1:t, L_vector)
plot(1:t, Lstar_vector,'--')
pbaspect([1 1 1])

legend('Present value of bank','Constrained optimal lending','Unconstrained optimal lending','Location','southeast')


x = [0.25 0.25];
y = [0.6 0.52];
p=annotation('textarrow',x,y,'String','Lending increase happens sooner ')
p.FontSize = 7;
xlabel('Time Period')
set(gca,'FontSize',7)
title('Short low interest rate period')
grid on
hold off
saveas(f3, strcat(data_dir, '\multi.eps'),'epsc')

