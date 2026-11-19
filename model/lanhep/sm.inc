% Standard Model of micrOMEGAs (STFM/work/lanhep/sm.inc) with a U(1)_X gauge boson
% kinetically mixed with hypercharge (Zp, exact Z-Zp mixing). Unitary gauge only.


external_func(alphaQCD,1).
external_func(McEff,1).
external_func(MbEff,1).
external_func(MtEff,1).
external_func(MqEff,2). 
external_func(initQCD5,4).
external_func(MixMatrix,3).
external_func(rDiagonal,*).
external_func(MassArray,2).

option ReduceGamma5=0.
let g5=gamma5.
use sm_tex.

parameter  EE  = 0.31333 : 'Electromagnetic coupling constant (<->1/128)',
	   GG  = 1.117   : 'Strong coupling constant (Z point)  (PDG-94)',
	   SW  = 0.4740  : 'sin of the Weinberg angle 0.474 - "on-shell",481 - "MS-bar" )'.


parameter  CW  = sqrt(1-SW**2) : 'cos of the Weinberg angle'.

% Masses 

parameter Q=100: 'Scale of effective running masses',
            MW=80.385: 'W boson mass',
            GF=EE**2/(2*SW*MW)**2/Sqrt2 :' experimental value 1.166E-5 [1/GeV^2]',  
            MZ0=MW/CW: 'Z mass without Z-Zp mixing',
            Mtp  =172.5:'Top quark pole mass',
            McMc =1.23:     'Mc(Mc)  MS-BAR',
            MbMb =4.25:     'Mb(Mb)  MS-BAR',
            alphaSMZ=0.1184:'Srtong alpha(MZ)', 
            LamQCD=initQCD5(alphaSMZ,McMc,MbMb,Mtp),
            Mb=MbEff(Q),
            Mc=McEff(Q),
            Ms=MqEff(0.096,Q): 's-quark effective mass  via 2MeV running one' .

%================== U(1)_X: neutral gauge boson masses and Z-Zp mixing
%  Inputs are the physical Zp mass MZp, the kinetic mixing epsX = sin(eps')
%  of -sin(eps')/2*B_{mu nu}X^{mu nu} and gX. The Z mass MZ and the U(1)_X
%  mass MX follow from the exact diagonalization of the (Zt,Ztp) mass matrix
%    {{MZ0**2, bZZ}, {bZZ, MZ0**2*SW**2*teps**2 + MX**2/ceps**2}},
%  with bZZ = MZ0**2*SW*teps and (Zt,Ztp) = {{cchi,-schi},{schi,cchi}}*(Z,Zp).
%  Valid for MZp below or above MZ0 = MW/CW, except in a narrow window just
%  above MZ0 where MX**2 < 0 (no solution); the window grows with epsX, e.g.
%  91.2921-91.2932 GeV for epsX = 1E-2. main.c stops there.

parameter  gX   = 0.01 : 'U(1)_X gauge coupling',
           epsX = 1E-6 : 'kinetic mixing, sin of eps prime',
           MZp  = 30   : 'Zp (dark photon) mass'.

parameter  ceps = sqrt(1-epsX**2) : 'cos of eps prime',
           teps = epsX/ceps       : 'tan of eps prime',
           bZZ  = MZ0**2*SW*teps : 'Z-Zp off-diagonal mass squared',
           dZZ  = MZp**2+bZZ**2/(MZ0**2-MZp**2) : 'Zp-Zp entry of mass squared matrix',
           MZ   = sqrt(MZ0**2+dZZ-MZp**2) : 'Z boson mass',
           MX   = ceps*sqrt(dZZ-MZ0**2*SW**2*teps**2) : '* U(1)_X mass before mixing',
           tchi = bZZ/(MZ0**2-MZp**2) : 'tan of Z-Zp mixing angle chi',
           cchi = 1/sqrt(1+tchi**2)   : 'cos of chi',
           schi = tchi*cchi           : '* sin of chi',
           epsDP= epsX*CW : '* dark photon coupling eps_DP = sin(eps prime)*CW'.


%  The Z-Zp mixing is written for the unitary gauge only.
do_if gauge_fixing==Feynman.
    write('Error: this model is written in the unitary gauge; set gauge_fixing=unitary.').
    quit.
end_if.

vector
	A/A: (photon, gauge),
	Z/Z:('Z boson', mass MZ, width wZ = auto),
	Zp/Zp:('Zp dark photon', pdg 32, mass MZp, width wZp = auto),
	G/G: (gluon, color c8, gauge),
        'W+'/'W-': ('W boson', mass MW, width wW = auto).

spinor 		ne/Ne:(neutrino,left),        e:(electron),
		nm/Nm:('mu-neutrino',left),   m:(muon),
		nl/Nl:('tau-neutrino',left),  l:('tau-lepton', mass Ml  = 1.777).

spinor		u:('u-quark',color c3, mass Mq = 0.05),
		d:('d-quark',color c3, mass Mq),
		c:('c-quark',color c3, mass Mc),
		s:('s-quark',color c3, mass Ms),
		t:('t-quark',color c3, mass Mtp, width wt = auto),
		b:('b-quark',color c3, mass Mb).


do_if CKMdim==1.
 let  Vub=0, Vcb=0, Vtd=0, Vts=0, Vtb=1, Vud=1, Vus=0, Vcs=1, Vcd=0.

do_else.

parameter  s12 = 0.221   : 'Parameter of C-K-M matrix (PDG-94)',
           c12  = sqrt(1-s12**2).

do_if CKMdim==2.

  parameter  Vud = c12  ,
             Vus = s12  ,
             Vcs = Vud  ,
             Vcd = -Vus .
  let        Vub = 0, Vcb = 0, Vtd = 0, Vts = 0, Vtb = 1.
  OrthMatrix({{Vud,Vus}, {Vcd,Vcs}}).

do_else.
   
  parameter  s23 = 0.040   : 'Parameter of C-K-M matrix (PDG-94)',
             s13 = 0.0035  : 'Parameter of C-K-M matrix (PDG-94)',
             c23  = sqrt(1-s23**2),
             c13  = sqrt(1-s13**2).
		  
parameter  Vud = c12*c13                ,
           Vus = s12*c13                ,
           Vub = s13     	        ,
           Vcd = (-s12*c23-c12*s23*s13) ,
           Vcs = (c12*c23-s12*s23*s13)  ,
           Vcb = s23*c13                ,
           Vtd = (s12*s23-c12*c23*s13) 	,
           Vts = (-c12*s23-s12*c23*s13)	,
           Vtb = c23*c13  	        .

OrthMatrix( { {Vud,Vus,Vub}, {Vcd,Vcs,Vcb}, {Vtd,Vts,Vtb}} ).

end_if.
end_if.




scalar h/h:(Higgs, mass Mh = 125, width wh =auto).

%================== dark Higgs and h-hX mixing (Higgs portal, see u1x.inc)
%  Gauge eigenstates: h0 in the doublet, s0 in phi_q = vq + s0/Sqrt2.
%  Mass eigenstates h (Mh) and hX (MhX):  h0 = ca*h - sa*hX,  s0 = sa*h + ca*hX.
%  sa = 0 is the model of the paper (no portal coupling).
parameter  MhX = 70  : 'dark Higgs mass',
           sa  = 0   : 'sin of h-hX mixing angle'.
parameter  ca  = sqrt(1-sa**2) : 'cos of h-hX mixing angle'.

scalar hX/hX:('dark Higgs', pdg 9000025, mass MhX, width whX = auto).

let h0 = ca*h - sa*hX.
									
let l1={ne,e}, L1={Ne,E}.
let l2={nm,m}, L2={Nm,M}.
let l3={nl,l}, L3={Nl,L}.

let q1={u,d}, Q1={U,D}, q1a={u,Vud*d+Vus*s+Vub*b}, Q1a={U,Vud*D+Vus*S+Vub*B}.
let q2={c,s}, Q2={C,S}, q2a={c,Vcd*d+Vcs*s+Vcb*b}, Q2a={C,Vcd*D+Vcs*S+Vcb*B}. 
let q3={t,b}, Q3={T,B}, q3a={t,Vtd*d+Vts*s+Vtb*b}, Q3a={T,Vtd*D+Vts*S+Vtb*B}.


%================== U(1)_Y and U(1)_X gauge fields in the mass basis (Z, Zp)
%  Zt, Ztp: neutral fields after removing the kinetic mixing, before Z-Zp mixing.
%  B1 is the hypercharge field (B tilde of the paper), XX the U(1)_X field.
let Zt  = cchi*Z - schi*Zp, Ztp = schi*Z + cchi*Zp.
let B1  = CW*A - SW*Zt - teps*Ztp.
let XX  = Ztp/ceps.

%================== SU(2) gauge field  and Faddeev-Popov ghosts 
let  W = { ('W+'+'W-')/Sqrt2,     i*('W+'-'W-')/Sqrt2,     CW*Zt+SW*A        }.
let  WW= { 'W-', CW*Zt+SW*A, 'W+'}.
do_if gauge_fixing==Feynman.
let  gh= { ('W+.c'+'W-.c')/Sqrt2, i*('W+.c'-'W-.c')/Sqrt2, CW*'Z.c'+SW*'A.c'}.
let  Gh= { ('W+.C'+'W-.C')/Sqrt2, i*('W+.C'-'W-.C')/Sqrt2, CW*'Z.C'+SW*'A.C'}. 
end_if.

%================== SU(2) and U(1) gauge couplings and VEV
let g2=EE/SW, g1=EE/CW.
let vevh=2*MW/EE*SW.      


%================== Self-interaction of gauge bosons

% Abelian kinetic terms with the kinetic mixing -sin(eps')/2*B_{mu nu}*X^{mu nu}, Eq. (13).
% They are diagonal in (A, Z, Zp) and give no vertices; kept for CheckMasses.
lterm -F**2/4  where  F=deriv^mu*B1^nu-deriv^nu*B1^mu.
lterm -F**2/4  where  F=deriv^mu*XX^nu-deriv^nu*XX^mu.
lterm -epsX/2*(deriv^mu*B1^nu-deriv^nu*B1^mu)*(deriv^mu*XX^nu-deriv^nu*XX^mu).
% photon gauge fixing (Feynman gauge for the massless photon, as CalcHEP assumes)
lterm -1/2*(deriv*A)**2.
lterm -F**2/4  where  F=deriv^mu*G^nu^a-deriv^nu*G^mu^a+i*GG*f_SU3^a^b^c*G^mu^b*G^nu^c.
lterm -F**2/4  where  F=deriv^mu*W^nu^a-deriv^nu*W^mu^a+g2*eps^a^b^c*W^mu^b*W^nu^c.


%================== left fermion interaction with gauge fields

lterm  	i*anti(psi)*gamma*(1-g5)/2*(deriv-i*g2*tau*W/2-i*Y*g1*B1)*psi
		where 
			psi=l1,  Y=-1/2;
			psi=l2,  Y=-1/2;
			psi=l3,  Y=-1/2;
			psi=q1a, Y= 1/6;
			psi=q2a, Y= 1/6;
			psi=q3a, Y= 1/6.

%================== right fermion interaction with gauge fields

lterm  	i*anti(psi)*gamma*(1+g5)/2*(deriv - i*Y*g1*B1)*psi
		where 
			psi=e,Y= -1;
			psi=m,Y= -1;
			psi=l,Y= -1;
			psi=u, Y=  2/3;
			psi=c, Y=  2/3;
			psi=t, Y=  2/3;
			psi=d, Y= -1/3;
			psi=s, Y= -1/3;
			psi=b, Y= -1/3.

%================== quark-gluon interaction

lterm  GG*anti(psi)*lambda*gamma*G*psi where
	psi=q1; psi=q2; psi=q3.


%=================  Higgs doublet
do_if gauge_fixing==Feynman.
    let shd = {  i*'W+.f',  (vev(vevh)+h-i*'Z.f')/Sqrt2 }.   % scalar Higgs doublet
do_else.
    let shd = { 0,          (vev(vevh)+h0)/Sqrt2        }. 
end_if.

let  shD = anti(shd).                                        % scalar Higgs anti-doublet. 

%===  Generalized dirivative of Higgs doublet 

let Dshd^mu^a = (deriv^mu-i*g1/2*B1^mu)*shd^a -i*g2/2*tau^a^b^c*W^mu^c*shd^b.
let DshD^mu^a = (deriv^mu+i*g1/2*B1^mu)*shD^a +i*g2/2*tau^b^a^c*W^mu^c*shD^b.
    
%===  Lagrangian 
lterm DshD*Dshd.     % the scalar potential, with the portal term, is in u1x.inc


%================ Yukawa terms  
%=== for bottom quarks 
lterm  -M*Sqrt2/vevh*anti(pl)*(1+g5)/2*pr*shd +  AddHermConj 
    where
	M=Vud*Mq, pl=q1a, pr=d;
	M=Vus*Ms, pl=q1a, pr=s;
	M=Vub*Mb, pl=q1a, pr=b;
	M=Vcd*0,  pl=q2a, pr=d;
	M=Vcs*Ms, pl=q2a, pr=s;
	M=Vcb*Mb, pl=q2a, pr=b;
	M=Vtd*0,  pl=q3a, pr=d;
	M=Vts*Ms, pl=q3a, pr=s;
	M=Vtb*Mb, pl=q3a, pr=b.

%===  for up quarks
   let eps2 = {{0,-1},{1,0}}.  % eps2^{12} =1 LanHEP reads matrices by columns 
   lterm   -M*Sqrt2/vevh*anti(pl)*(1+g5)/2*eps2*pr*shD +  AddHermConj
   where
	M=Mq,  pl=q1a, pr=u;
	M=Mc,  pl=q2a, pr=c;
	M=Mtp, pl=q3a, pr=t.

%=== for charged leptons 
lterm  -M*Sqrt2/vevh*anti(pl)*(1+g5)/2*pr*shd +  AddHermConj 
    where
%	M=Mm,   pl=l2,  pr=m;
	M=Ml,  pl=l3,  pr=l.
	
%====================== gauge fixing terms and corresponding ghost terms

%=== for SU(3)
lterm  -1/2*(deriv*G)**2.
lterm -i*GG*f_SU3*ccghost(G)*G^mu*deriv^mu*ghost(G).


%=== for SU(2)*U(1) 

do_if gauge_fixing==Feynman.

lterm  -1/2*(deriv*A)**2  -1/2*(2*(deriv*'W+'+vevh*g2/2*'W+.f')*(deriv*'W-'+vevh*g2/2*'W-.f') +(deriv*Z+ vevh*g2/2/CW*'Z.f')**2).
                
lterm -g2*eps*deriv*Gh*gh*W.

lterm vevh*g2**2/4*(-h*('W-.C'*'W+.c' + 'W+.C'*'W-.c' + 'Z.C'*'Z.c'/CW**2)  +i*'Z.f'*('W-.C'*'W+.c' - 'W+.C'*'W-.c') ).

lterm -i*g2**2*vevh/4/CW*(
	'W+.f'*('W-.C'*'Z.c'*(1-2*SW**2)+'W-.c'*'Z.C' + 2*CW*SW*'W-.C'*'A.c')
      -	'W-.f'*('W+.C'*'Z.c'*(1-2*SW**2)+'W+.c'*'Z.C' + 2*CW*SW*'W+.C'*'A.c')  
                    ).
end_if.


SetAngle(1-SW**2=CW**2).
SetEM(A,EE).


