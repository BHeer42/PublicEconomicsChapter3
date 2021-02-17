@ ----------------------------  Ch3_olg_dyn1.g -------------------------------------

computes the dynamics of the Example economy

author: Burkhard Heer

this version: Feb 13, 2013

------------------------------------------------------------------------------ @

new; clear all; cls;
library pgraph; //, user;
MachEps=machEpsilon;
#include Derivatives.src;

_plotstyle=1;     // 0 -- figure generated with xy(x,y), GAUSS 13,.. , 1 -- figure generated with plotXY(.,.,.) using GAUSS 19,...

alpha=0.36;		// production elasticity of capital

beta1=0.40;		// discount factor 

n=0.1;			// population growth rate

k=( beta1/(1+beta1) * (1-alpha)/(1+n) )^(1/(1-alpha));	// steady state value of capital stock 

k0=k/3;		// initial capital stock 

tt=seqa(0,1,20);	// periods

kt=zeros(20,1);		// time series for capital stock

kt[1]=k0;
i=1;
do until i==20;
	i=i+1;
	kt[i]=beta1/(1+beta1) * (1-alpha)/(1+n) * kt[i-1]^(alpha);
endo;


if _plotstyle==0;
    set_xy;
    //title("Dynamics of the capital stock");
    xlabel("Period t");
    yfont= "Capital stock k\201]t[";
    ylabel(yfont);
    xy(tt,kt);
else;
    struct PlotControl myPlot;
	myplot=plotGetDefaults("xy");	
    plotSetXLabel(&myplot,"Period t");
    plotSetYLabel(&myplot,"Capital Stock");
    plotSetTitle(&myPlot,"Dynamics of the Capital Stock");
    plotXY(myplot,tt,kt);
endif;                

wait;
"comutation of numerical derivative:";
Jac=CDJac(&kdyn,k,1);
Jac;

@ ---------------------------------------------  procedures --------------------------------- @

proc kdyn(x);
	retp( beta1/(1+beta1) * (1-alpha)/(1+n) * x^(alpha) );
endp;


proc(0)=set_xy();

    graphset;
    _pmcolor = 0|0|0|0|0|0|0|0|15; 
    _pframe = { 1, 1};
    _pcolor = { 12, 2, 1, 0}; 
    _pltype = { 6, 3, 5, 6 };
    _pcross=0;
	_pgrid = { 2, 2}; 
    _ptitlht = 0.25;
    _plwidth = 10;
	_plwidth=8;
    _paxht=0.25;
    _pnumht=0.20;
    _pdate=0;
    fonts("Microb");
    fonts("simplex complex microb simgrma");
    pause(0.5);

retp();
endp;
