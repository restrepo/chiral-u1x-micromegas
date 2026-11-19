% Chiral U(1)_X dark sector of arXiv:2506.04718 (He, Nomura, Yokozaki),
% benchmark (Q1,Q2,QM) = (-9,-1,-5), dark Higgs charge q = 10.
%
% Five left-handed Weyl fermions with X charges -9, -5, -1, 7, 8.
%   ~X  : Dirac DM,      P_L ~X = psi_{-1,L},  P_R ~X = (psi_{-9,L})^c  [Z2 odd]
%   ~~M : Majorana DM,   P_L ~~M = psi_{-5,L}                         [Z2' odd]
%   nX  : massless Dirac (dark radiation), P_L nX = psi_{7,L}, P_R nX = (psi_{8,L})^c
% The two massless Weyl fermions psi_7 and psi_8 are written as one massless
% Dirac field; with no mass term this is exact.
% Unitary gauge: H = (0, (v + h0)/Sqrt2), phi_q = vq + s0/Sqrt2, with
%   V = lamH*(|H|**2 - v**2/2)**2 + lamX*(|phi_q|**2 - vq**2)**2
%     + lamHX*(|H|**2 - v**2/2)*(|phi_q|**2 - vq**2).
% The paper has lamHX = 0. Here the inputs are Mh, MhX and sa (h-hX mixing),
% so MX**2 = 2*gX**2*q**2*vq**2, MD = yD*vq, MM = yM*vq.

let qphi = 10.                     % U(1)_X charge of the dark Higgs, q = -(Q1+Q2) = -2*QM
let QD1 = -9, QD2 = -1, QMj = -5.  % Dirac pair and Majorana charge
let QR1 = 7,  QR2 = 8.             % massless fermions

parameter  MDm = 14.9  : 'Dirac DM mass MD = yD*vq',
           MMj = 40    : 'Majorana DM mass MM = yM*vq'.

parameter  vq   = MX/(Sqrt2*gX*10) : '* dark Higgs VEV, MX**2 = 2*gX**2*q**2*vq**2',
           lamH = (ca**2*Mh**2+sa**2*MhX**2)/(2*(2*MW*SW/EE)**2) : 'Higgs quartic coupling',
           lamX = (sa**2*Mh**2+ca**2*MhX**2)/(4*vq**2)  : 'dark Higgs quartic coupling',
           lamHX= ca*sa*(Mh**2-MhX**2)/(Sqrt2*(2*MW*SW/EE)*vq) : '* Higgs portal coupling'.

spinor  '~X'/'~x':('Dirac DM',    pdg 9000011, mass MDm),
        '~~M'/'~~M':('Majorana DM', pdg 9000012, mass MMj),
        nX/NX:('dark radiation', pdg 9000013).

%================== gauge interactions of the dark fermions
%  D_mu = deriv_mu - i*gX*Q*X_mu for a field of charge Q; the right-handed
%  component of a Dirac field built from (psi_{Q,L})^c has charge -Q.

lterm  i*anti(psi)*gamma*( (1-g5)/2*(deriv - i*gX*QL*XX)
                          +(1+g5)/2*(deriv - i*gX*QR*XX) )*psi
   where
     psi='~X', QL=QD2, QR=-QD1;
     psi=nX,   QL=QR1, QR=-QR2.

%  Majorana: (1/2)*anti(psiM)*i*gamma*(deriv - i*gX*XX*(-QM)*g5)*psiM, Eq. (10)
lterm  i/2*cc('~~M')*gamma*deriv*'~~M'
      + gX/2*(-QMj)*cc('~~M')*gamma*g5*XX*'~~M'.

%================== Yukawa terms: masses and dark Higgs couplings
%  yD*phi_q*psi_{Q1}^c psi_{Q2} + yM/2*phi_q*psi_{QM}^c psi_{QM} + h.c.
let s0 = sa*h + ca*hX.

lterm  -MDm*(1+s0/Sqrt2/vq)*anti('~X')*'~X'.
lterm  -MMj/2*(1+s0/Sqrt2/vq)*cc('~~M')*'~~M'.

%================== dark Higgs
let phiX = vev(vq) + s0/Sqrt2.

lterm  deriv^mu*phiX*deriv^mu*phiX
      + gX**2*qphi**2*XX^mu*XX^mu*phiX**2.

%================== scalar potential with the Higgs portal
lterm  - lamH*(shD*shd - vevh**2/2)**2
       - lamX*(phiX**2 - vq**2)**2
       - lamHX*(shD*shd - vevh**2/2)*(phiX**2 - vq**2).
