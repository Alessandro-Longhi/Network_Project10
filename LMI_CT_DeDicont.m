function [K,rho,feas]=LMI_CT_DeDicont(A,B,C,N,ContStruc,type,value)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - A: system matrix.
% - B: input matrices (i.e., B{1},..., B{N} are the input matrices of the decomposed system, one for each channel).
% - C: output matrices  (i.e., C{1},..., C{N} are the output matrices of the decomposed system, one for each channel, where [Cdec{1}',...,
% Cdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral abscissa of matrix (A+B*K) - note that [C{1}',...,
% C{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)

Btot=[];
for i=1:N
    m(i)=size(B{i},2);
    n(i)=size(C{i},1);
    Btot=[Btot,B{i}];
end
ntot=size(A,1);
mtot=sum(m);

yalmip clear

if ContStruc==ones(N,N)
    % Centralized design
    Y=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Decentralized/distributed design
    Y=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        Y=blkdiag(Y,sdpvar(n(i)));
        ninc=0;
        for j=1:N
            if ContStruc(i,j)==0
                L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
            end
            ninc=ninc+n(j);
        end
        minc=minc+m(i);
    end  
end

LMIconstr=[Y*A'+A*Y+Btot*L+L'*Btot'<=-1e-2*eye(ntot)]+[Y>=1e-2*eye(ntot)];

% LEGEND
% 1 - LMI for performance
%
if type == 1
    a = value;
    LMIconstr = LMIconstr + [Y*A'+A*Y+L'*Btot'+Btot*L+2*a*Y<=-1e-2*eye(ntot)];
    objective = [];
end
if type == 2
    theta = value;
    LMIconstr = LMIconstr + [[sind(theta)*(A*Y+Y*A'+Btot*L+L'*Btot') cosd(theta)*(A*Y-Y*A'+Btot*L-L'*Btot');
                              cosd(theta)*(-A*Y+Y*A'-Btot*L+L'*Btot') sind(theta)*(A*Y+Y*A'+Btot*L+L'*Btot')]<=-1e-2*eye(ntot*2)];
    objective = [];
end

if type == 3
    a = value;
    aL = 0.01;
    aY = 10;
    KY=sdpvar(1);
    KL=sdpvar(1);
    I = eye(ntot);
    objective = aY*KY + aL*KL;
    LMIconstr = LMIconstr + [[KL*eye(ntot) L'; L eye(mtot)]>=1e-2*eye(ntot+mtot)] + [[KY*I I; I Y]>=1e-2*eye(ntot*2)]+ [Y*A'+A*Y+L'*Btot'+Btot*L+2*a*Y<=-1e-2*eye(ntot)];
end

if type == 4
    Hw = [eye(19);
            zeros(2,19)];
    Dw = [zeros(19,2);
            eye(2,2)];
    Bw = eye(19);
    S=sdpvar(21);
    objective = trace(S);
    LMIconstr = LMIconstr + [[Y*A'+A*Y+Btot*L+L'*Btot'+Bw*Bw'] <= -1e-2*eye(ntot)] + [[S Hw*Y+Dw*L;
                                                                                 L'*Dw'+Y*Hw' Y]>= 1e-2*eye(ntot*2+mtot)];
end

options=sdpsettings('solver','sedumi');
J=optimize(LMIconstr,objective,options);
feas=J.problem;
L=double(L);
Y=double(Y);
if type == 3
    KL=double(KL)
    KY=double(KY)
end
K=L/Y;
rho=max(real(eig(A+Btot*K)));
