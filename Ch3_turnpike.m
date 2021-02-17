% ----------------------------  Ch3_turnpike.m ----------------------------
%
% computes the dynamics of the central planner OLG economy 
% from Chapter 3.3.4 in my Book
%
% author: Burkhard Heer
%
% this version: June 6, 2018
%
% -------------------------------------------------------------------------



clear; clc;
%parameters


par.alpha=0.36;		% production elasticity of capital
par.beta=0.4;			%  discount factor 

par.n=0.1;			% population growth rate

k=(par.alpha/par.n)^(1/(1-par.alpha));	% steady state value of capital stock 

par.k0=0.7*k;		% initial capital stock 

par.kfinal=0.7*k;   % final capital stock

par.bigt=20;        % number of transition periods

tt=0:1:par.bigt+2;	% periods

kt=zeros(par.bigt+2,1);		% time series for capital stock

x=kdyn(par,k);

kguess=1.1*par.k0;          % initial guess for solution

ksol = fsolve(@(k)kdyn(par,k),kguess);


    kt2=zeros(par.bigt+2,1);
    kt2(1)=par.k0;
	kt2(2)=ksol;
	for i=3:par.bigt+2,
		kt2(i)=kt2(i-2)+kt2(i-2)^(par.alpha) -(1+par.n)*kt2(i-1) ;
        kt2(i)=(kt2(i-1)+kt2(i-1)^(par.alpha))/(1+par.n) - (1+par.alpha*kt2(i-1)^(par.alpha-1))/((1+par.n)^2) * kt2(i);
    end;


figure
plot([0:par.bigt+1],kt2); hold on
xlabel('Period t')
ylabel('Capital stock k_t')
