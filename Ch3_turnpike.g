@ ----------------------------  Ch3_turnpike.g -------------------------------------

computes the dynamics of the central planner OLG economy from Chapter 1 in my OLG class

also demonstrates that Alfred Maußner's non-linear equation solving routine "FixvMN1"
is able to find the solution, while the routine "eqSolve" is unable to find the solution!!

author: Burkhard Heer

this version: Feb 15, 2013

------------------------------------------------------------------------------ @


new; clear all; cls;
library pgraph, user;

MachEps=machEpsilon;
#include Derivatives.src;
#include toolbox.src;


_plotstyle=1;     // 0 -- figure generated with xy(x,y), GAUSS 13,.. , 1 -- figure generated with plotXY(.,.,.) using GAUSS 19,...


alpha=0.36;		// production elasticity of capital
b=0.4;			//  discount factor 

n=0.1;			// population growth rate

k=(alpha/n)^(1/(1-alpha));	// steady state value of capital stock 

beta1=0.40;		//
k0=0.7*k;		// initial capital stock 

kfinal=0.7*k;

bigt=20;
//bigt=100;

tt=seqa(0,1,bigt+2);	// periods

kt=zeros(bigt+2,1);		// time series for capital stock


ygrid=zeros(1000,1);
kgrid=seqa(k0,(k-k0)/1001,1000);
i=0; do until i==1000; i=i+1; ygrid[i]=kdyn(kgrid[i]); endo;
temp=kgrid~ygrid;
temp1=delif(temp,temp[.,2].==miss(1,1));


set_xy;
xlabel("Period t");
yfont= "kdyn";
ylabel(yfont);
xy(temp1[.,1],temp1[.,2]);
wait;


k1=k0;		// initial guess for k_1


{k1,crit}=FixVMN1(k1,&kdyn);
if crit[1]/=0; "no convergence of FixVMN1"; wait; endif;

"retcode FixVMN1: " crit;
wait;
"kdyn(k1)";
kdyn(k1);
wait;

if crit[1]==0;

        
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
endif;

k1=1.1*k0;		// initial guess for k_1
 {k1,crit}=eqSolve(&kdyn,k1);
if crit/=1; "no convergence of eqSolve"; wait; endif;

"retcode eqSolve: " crit;
wait;
"kdyn(k1)";
kdyn(k1);
wait;


if crit==0;
	
	
set_xy;
//title("Dynamics of the capital stock");
xlabel("Period t");
yfont= "Capital stock k\201]t[";
ylabel(yfont);
xy(tt,kt);
	
endif;

proc kdyn(x);
	local i;
	
	kt[1]=k0;
	kt[2]=x;
	i=2;
	do until i==bigt+2;
		i=i+1;
		kt[i]=(kt[i-1]+kt[i-1]^(alpha) ) / (1+n) - (1+alpha*kt[i-1]^(alpha-1)) / (1+n)^2 * ( kt[i-2]+kt[i-2]^(alpha) -(1+n)*kt[i-1] );
	endo;
	kt[bigt+2]-kfinal; 
	retp(kt[bigt+2]-kfinal);
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
