#/usr/bin/env python
from __future__ import division
import numpy as np
global alpha,mc,mb,mt,mw,mz,mu,muw,alphasz,z,lambda1,lambda2
import scipy.special as special
from mpmath import polylog as PolyLog
from mpmath import mp as mp
from scipy.integrate import quad as Integrate
mp.dps=15;mp.pretty=True

alpha=1/137.056;
mc=1.27;
mb=4.8;
mt=164;
mw=80.399;
mz=91.1876;
mu=mb;
muw=mw;
alphasz=0.119;
z=0.31;
#mathematica compatibility functions:
Log=np.log
Pi=np.pi
I=1j
Zeta=lambda x: special.zetac(x)+1
ArcTan=np.arctan
Sqrt=np.sqrt
Abs=np.abs

#funcion cinematica 
fz=1-8*z**2+8*z**6-z**8-24*z**4*Log(z);
#(* strong coupling constant *)
v=lambda y: 1-23/3*alphasz/(2*Pi)*Log(mz/y);
alphas=lambda y:alphasz/v(y)*(1-116/23*alphasz/(4*Pi)*Log(v(y))/v(y));
#(* operators at mw scale *)
xtw=(mt/mw)**2;
yth=lambda mH: (mt/mH)**2;

#(* mH=500; U2=1; U1=1; *)
c20=1;
c7sm=lambda x:(3*x**3-2*x**2)/(4*(x-1)**4)*Log(x)+(-8*x**3-5*x**2+7*x)/(24*(x-1)**3);
c8sm=lambda x:(-3*x**2)/(4*(x-1)**4)*Log(x)+(-x**3+5*x**2+2*x)/(8*(x-1)**3);
c7yy=lambda x:x*((-8*x**2-5*x+7)/(72*(x-1)**3)+(3*x**2-2*x)/(12*(x-1)**4)*Log(x));
c7xy=lambda x:x/12*((-5*x**2+8*x-3+(6*x-4)*Log(x))/(x-1)**3);
c8yy=lambda x:( -x**3+5*x**2+2*x)/(24*(x-1)**3)-(3*x**2)/(12*(x-1)**4)*Log(x);
c8xy=lambda x:x/4*((- x**2+4*x-3-2* Log(x))/(x-1)**3);
c70=lambda x,y,U1,U2:c7sm(x)+U1**2*c7yy(y)+U2*U1*c7xy(y);
c80=lambda x,y,U1,U2:c8sm(x)+U1**2*c8yy(y)+U2*U1*c8xy(y); 

#(* effective operator c7 *)
eta=alphas(mw)/alphas(mu);
c7ef=lambda x,y,U1,U2:eta**(16/23)*c70(x,y,U1,U2)+8/3*(eta**(14/23)-eta**(16/23))*c80(x,y,U1,U2)\
      +c20*(2.2996*eta**(14/23)-1.0880*eta**(16/23)-3/7*eta**(6/23)\
        -1/14*eta**(-12/23)-0.6494*eta**0.4086-0.0380*eta**-0.4230-0.0185*eta**\
        -0.8994-0.0057*eta**0.1456);
#(* the branching is by def *)
#(* ((Vts*Vtb)/Vcb)**2 is CKM factor included numerically and the semileptonic ratio from pdg *)
BrLO=lambda mH,U1,U2:10.08/100*0.971*(6*alpha)/(Pi*fz)*(c7ef(xtw,yth(mH),U1,U2))**2;
print BrLO(150, 0, 0)


#(* NLO *)
lambda1 =0.5;
lambda2=-0.12;
kz=1-(2*alphas(mu))/(3*Pi)*((Pi**2-31/4)*(1-z)**2+1.5)+1/mb**2*(lambda1/2+(3*lambda2)/2*(1-4*(1-z**2)**4/fz));
c1ef=eta**(6/23)-eta**(-12/23);
c2ef=2/3*eta**(6/23)+1/3*eta**(-12/23);
c8ef=lambda mH,U1,U2:eta**(14/23)*c80(xtw,yth(mH),U1,U2)+0.8623*eta**(14/23)-0.9135*eta**0.4086+0.0873*eta**-0.4230-0.0571*eta**-0.8994+0.0209*eta**0.1456;
fg=0.29**2;
r2=2/243*(-833+144*Pi**2*fg**1.5+(1728-180*Pi**2-1296*Zeta(3)+(1296-324*Pi**2)*Log(fg)+108*Log(fg)**2+36*Log(fg)**3)*fg+(648+72*Pi**2+(432-216*Pi**2)*Log(fg)+36*Log(fg)**3)*fg**2+(-54-84*Pi**2+1092*Log(fg)-756*Log(fg)**2)*fg**3)+(16*Pi*I)/81*(-5+(45-3*Pi**2+9*Log(fg)+9*Log(fg)**2)*fg+(-3*Pi**2+9*Log(fg)**2)*fg**2+(28-12*Log(fg))*fg**3);
r1=-(1/6)*r2;
r7= 32/9-8/9*Pi**2;
r8=-(4/27)*(-33+2*Pi**2-6*Pi*I);
vmu= lambda mH,x,y,U1,U2: alphas(mu)/(4*Pi)*(c1ef*r1+c2ef*r2+c7ef(xtw,yth(mH),U1,U2)*r7+c8ef(mH,U1,U2)*r8 -16/3*c7ef(xtw,yth(mH),U1,U2));
e0=lambda x: (x*(x**2+11*x-18))/(12*(x-1)**3)+(x**2*(4*x**2-16*x+15)*Log(x))/(6*(x-1)**4)-2/3*Log(x)-2/3;
eh=lambda x:x/36*((7*x**3-36*x**2+45*x-16 +(18*x-12)*Log(x))/(x-1)**4);
c41=lambda mH,U1,U2:e0(xtw)+ 2/3*Log(muw**2/mw**2)+U1**2*eh(yth(mH));
w7sm=lambda x:(-16*x**4-122*x**3+80*x**2-8*x)/(9*(x-1)**4)*PolyLog(2,1-1/x) +(6*x**4+46* x**3-28*x**2)/(3*(x-1)**5)*Log(x)**2+1/(81*(x-1)**5)*(-102*x**5-588*x**4-2262*x**3+3244*x**2-1364*x +208)*Log(x) +(1646* x**4+12205*x**3-10740*x**2+2509*x-436)/(486*(x-1)**4);
m7sm=lambda x:1/(81*(x-1)**5)*(82*x**5+301*x**4+703*x**3-2197*x**2+1319*x-208-(162*x**4+1242*x**3-756*x**2)*Log(x));
t7sm=lambda x:x/3*(1/(x-1)**5*(47*x**3-63*x**2+9*x+7-(18*x**3+30*x**2-24*x)*Log(x)) );
c71sm=lambda  x:w7sm(x)+m7sm(x)*Log(muw**2/mw**2)+t7sm(x)*(Log(mt**2/muw**2)-4/3);
w7yy=lambda x:(2*x)/9*((8*x**3-37*x**2+18*x)/(x-1)**4*PolyLog(2,1-1/x)+(3*x**3+23*x**2-14*x)/(x-1)**5*Log(x)**2+(21*x**4-192*x**3-174*x**2+251*x -50)/(9*(x-1)**5)*Log(x)+(-1202*x**3+7569*x**2-5436*x+797)/(108*(x-1)**4) ) -4/9*eh(x);
m7yy=lambda x:x/(27*(x-1)**5)*(-14*x**4+149*x**3-153*x**2-13*x+31-(18*x**3+138*x**2-84*x)*Log(x));
c71yy=lambda mH,x: w7yy(x)+m7yy(x)*Log(muw**2/mH**2)+t7sm(x)/3*(Log(mt**2/muw**2)-4/3);
w7xy=lambda x:4/3*((8*x**3-28*x**2+12*x)/(3*(x-1)**3)*PolyLog(2,1-1/x)+(3*x**3+14*x**2-8*x)/(3*(x-1)**4)*Log(x)**2+(4*x**4-24*x**3+2*x**2+6*x )/(3*(x-1)**4)*Log(x)+(-2*x**3+13*x**2-7*x)/((x-1)**3) ) ;
m7xy=lambda x:2/(9*(x-1)**4)*(-8*x**4+55*x**3-68*x**2+21*x-(6*x**3+28*x**2-16*x)*Log(x));
t7xy=lambda x:(2*x)/3*((13*x**2-20*x+7-(6*x**2+4*x-4)*Log(x))/(x-1)**4 );
c71xy=lambda mH,x:w7xy(x)+m7xy(x)*Log(muw**2/mH**2)+t7xy(x)*(Log(mt**2/muw**2)-4/3);

c71=lambda mH,U1,U2: c71sm(xtw)+U1**2*c71yy(mH,yth(mH))+U1*U2*c71xy(mH,yth(mH)) ;
w8sm=lambda x:(-4*x**4+40*x**3+41*x**2+x)/(6*(x-1)**4)*PolyLog(2,1-1/x)+(-17*x**3-31*x**2)/(2*(x-1)**5)*Log(x)**2+(-210*x**5+1086*x**4+4893*x**3+2857*x**2-1994*x+280)/(216*(x-1)**5)*Log(x)+(737*x**4-14102*x**3-28209*x**2+610*x-508)/(1296*(x-1)**4);
m8sm=lambda x: 1/(108*(x-1)**5)*(77*x**5-475*x**4-1111*x**3+607*x**2+1042*x-140+(918*x**3+1674*x**2)*Log(x));
t8sm=lambda x:2*x*((-x**3-9*x**2+9*x+1+(6*x**2+6*x)*Log(x))/(x-1)**5 ) ;
w8yy=lambda x: x /6*((13*x**3-17*x**2+30*x)/(x-1)**4*PolyLog(2,1-1/x)-(17*x**2+31*x)/(x-1)**5*Log(x)**2+(42*x**4+318*x**3+1353*x**2+817*x -226)/(36*(x-1)**5)*Log(x)+(-4451*x**3+7650*x**2-18153*x+1130)/(216*(x-1)**4) ) -1/6*eh(x) ;
m8yy=lambda x: x/(36*(x-1)**5)*(-7*x**4+25*x**3-279*x**2+223*x+38+(102*x**2+186*x)*Log(x));
w8xy=lambda x:1/3*((17*x**3-25*x**2+36*x)/(2*(x-1)**3)*PolyLog(2,1-1/x)-(17*x**2+19*x)/ ((x-1)**4)*Log(x)**2+(14*x**4-12*x**3+187*x**2+3*x )/(4*(x-1)**4)*Log(x)-(3*(29*x**3-44*x**2+143*x))/( 8*(x-1)**3) );
m8xy=lambda x:1/(6*(x-1)**4)*(-7*x**4+23*x**3-97*x**2+81*x+(34*x**2+38*x)*Log(x));
t8xy=lambda x: 2*x*((-x**2-4*x+5+(4*x+2)*Log(x))/(x-1)**4 );
c81sm=lambda x:w8sm(x)+m8sm(x)*Log(muw**2/mw**2)+t8sm(x)*(Log(mt**2/muw**2)-4/3) ;
c81yy=lambda mH,x:w8yy(x)+m8yy(x)*Log(muw**2/mH**2)+t8sm(x)/3*(Log(mt**2/muw**2)-4/3);
c81xy=lambda mH,x:w8xy(x)+m8xy(x)*Log(muw**2/mH**2)+t8xy(x)*(Log(mt**2/muw**2)-4/3);
c81=lambda mH,U1,U2:c81sm(xtw)+U1**2*c81yy(mH,yth(mH))+U1*U2*c81xy(mH,yth(mH));


c11ef=15+6*Log(muw**2/mw**2) ;
c71eff=lambda mH,U1,U2:eta**(39/23)*c71(mH,U1,U2)+8/3*(eta**(37/23)-eta**(39/23))*c81(mH,U1,U2)+c80(xtw,yth(mH),U1,U2)*(-(7164416/357075)*eta**(14/23)+297664/14283*eta**(16/23)+256868/14283*eta**(37/23)-6698884/357075*eta**(39/23))+37208/4761*(eta**(39/23)-eta**(16/23))*c70(xtw,yth(mH),U1,U2)+(4661194/816831*eta**(14/23)-8516/2217*eta**(16/23)-1.9043*eta**0.4086-0.1008*eta**-0.4230+0.1216*eta**-0.8994+0.0183*eta**0.1456)*eta*c41(mH,U1,U2) -17.3023*eta**(14/23)+8.5027*eta**(16/23)+4.5508*eta**(6/23)+0.7519*eta**(-12/23)+2.0040*eta**0.4086+0.7476*eta**-0.4230-0.5385*eta**-0.8994+0.0914*eta**0.1456+(9.9372*eta**(14/23)-7.4878*eta**(16/23)+1.2688*eta**(6/23)-0.2925*eta**(-12/23) -2.2923*eta**0.4086-0.1461*eta**-0.4230+0.1239*eta**-0.8994+0.0812*eta**0.1456 )*eta+(0.5784*eta**(14/23)-0.3921*eta**(16/23)-0.1429*eta**(6/23)+0.0476*eta**(-12/23) -0.1275*eta**0.4086+0.0317*eta**-0.4230+0.0078*eta**-0.8994-0.0031*eta**0.1456)*eta*c11ef;
c7eft=lambda mH,U1,U2:c7ef(xtw,yth(mH),U1,U2)+ alphas(mu)/(4*Pi)*c71eff(mH,U1,U2);
bd=lambda mH,U1,U2:c7eft(mH,U1,U2)+vmu(mH,xtw,yth(mH),U1,U2);

def Gt(t):
    if t<4:
        return -2*ArcTan(Sqrt(t/(4-t)))**2
    else:
        return -(Pi**2/2)+2*Log((Sqrt(t)+Sqrt(t-4))/2)**2-2*Pi*I*Log((Sqrt(t)+Sqrt(t-4))/2)

f220=(16*fg)/27*Integrate(lambda t: (1-fg*t)**2*Abs(Gt(t)/t+0.5)**2,0,3.999)[0];

f22p=(16*fg)/27*Integrate(lambda t: (1-fg*t)**2*Abs(Gt(t)/t+0.5)**2,4,11.89)[0];

f22=f220+f22p;

f270=-((8*fg**2)/9)*Integrate(lambda t: (1-fg*t)*(Gt(t)+0.5*t),0,3.999)[0]
f27p=-((8*fg**2)/9)*Integrate(lambda t: (1-fg*t)*(Gt(t)+0.5*t),4,11.89)[0]

f27=f270+f27p;
f11=f22/36;
f12=-(f22/3);
f17=-(f27/6);
f28=-(f27/3);
f18=-(f28/6);
f78=8/9*(25/12-Pi**2/6);
f88=1/27*(16/3-(4*Pi**2)/3);

A=lambda mH,U1,U2:alphas(mu)/ Pi*(c1ef**2*f11+c1ef*c2ef*f12+c1ef*c7ef(xtw,yth(mH),U1,U2)*f17+c1ef*c8ef(mH,U1,U2)*f18+c2ef**2*f22+c2ef*c7ef(xtw,yth(mH),U1,U2)*f27+c2ef*c8ef(mH,U1,U2)*f28+c7ef(xtw,yth(mH),U1,U2)*c8ef(mH,U1,U2)*f78+c8ef(mH,U1,U2)**2*f88);

Delta=lambda mH,U1,U2:1/mb**2*(lambda1*0.5-4.5*lambda2)*c7ef(xtw,yth(mH),U1,U2)**2+(-lambda2/(9*mc**2))*(c7ef(xtw,yth(mH),U1,U2)*(c2ef-c1ef/6));

brNLO=lambda mH,U1,U2:(10.08/100*0.971*(6*alpha)/(Pi*fz*kz )*(Abs(bd(mH,U1,U2))**2+ A(mH,U1,U2)+ Delta(mH,U1,U2))).real;

print brNLO(150,0,0)
other='''


N(brNLO(150,0,0))
'''


