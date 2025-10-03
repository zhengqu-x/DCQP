function [c,val,sdp_status,S] = generate_cut_dnn(Q,d,A,b,Aeq,beq,m,n,nuR,barx,x0,tol_mosek,beta)

% ======================================================================= %
% Generate cut for the following problem
%     min   x'Qx+2d'x
%     s.t.  Ax <= b
%           Aeq = beq
% at KKT point barx.
% 
% Solve the SDP problem
%     min c'*(x0-barx)
%     s.t. [Q d; d' -nuR]=S+[-A b; zeros(1,n) 1]'T[-A b; zeros(1,n) 1]
%                              +1/2[Q\bar x+d;-barx'*Q*barx-d'*barx+beta]*[-c' 1+c'*barx]
%                              +1/2[-c; 1+c'*barx]*[(Q\bar x+d)'
%                              -barx'*Q*bar-d'*barx+beta]+ U*[Aeq beq]+[Aeq beq]'*U'.

% INPUT
%
% Q, d             Parameters of the quadratic objective function
% A, b             Parameters of the inequality constraint A x <= b
% Aeq, beq         Equality constraints: Aeq x = beq
% m                Size(A,1)
% n                Size(A,2)
% nuR              \nu_R
% barx             \bar x
% x0               x_0
% tol_mosek        Tolerance parameter used in MOSEK
% beta             Beta
%
% OUTPUT
%
% c                c
% val              Primal optimal value returned by MOSEK
% sdp_status       Status of MOSEK
% S                S
% ======================================================================= %   



[~, res] = mosekopt('symbcon echo(0)');
symbcon = res.symbcon;

q=[Q*barx+d;-barx'*Q*barx-d'*barx+beta];
barQ=[Q d; d' -nuR];
barA=[-A b;zeros(1,n) 1];  % barA of dimension (m+1,n+1) 


% Solve the following SDP
%                                min  c'*(x0-barx)
%                                s.t. barQ-barA'*T*barA+0.5q*[c' -c'*barx-1]+0.5*[c;-1-c'*barx]*q'+ U*[Aeq beq]+[Aeq beq]'*U'=S
%                                     T>=0, S PSD.

% T_{1,1}=...=T_{m+1,m+1}=0
% Vec(T)=[T_{12};...;T_{1,m+1};T_{23};....;T_{2,m+1};....T_{m-1,m+1};T_{m,m+1}]
% sVec(barA'*T*barA)= M_T* Vec(T)

M_T=zeros((n+1)*(n+2)/2,(m+1)*m/2);   
idx=0;
for i=1:m
    for j=i+1:m+1
        idx=idx+1;
        Ai=barA(i,:);
        Aj=barA(j,:);
        Qij=Ai'*Aj;
        M_T(:,idx)=sVec(Qij+Qij');
    end
end
M_T=sparse(M_T);


% sVec(0.5q*[c' -c'*barx]+0.5*[c;-c'*barx]*q')=M_c*c
M_c=zeros((n+1)*(n+2)/2,n);
for i=1:n
    ei=zeros(n,1);
    ei(i)=1;
    tmp=q*[ei' -barx(i)];
    M_c(:,i)=0.5*sVec(tmp+tmp');
end

% Rq=sVec(0.5*q*[0' -1]+0.5*[0;-1]*q')

tmp=zeros(n+1,1);
tmp(end)=-1;
Rq=0.5*sVec(q*tmp'+tmp*q');


M_U=[];
meq=size(Aeq,1);
if meq>0
    zm=zeros(n+1,meq);
    idx=0;
    M_U=zeros((n+1)*(n+2)/2,(n+1)*meq);
    for i=1:n+1
        for j=1:meq
            idx=idx+1;
            Eij=zm;
            Eij(i,j)=1;
            M_U(:,idx)=sVec(Eij*[Aeq -beq]+[Aeq -beq]'*Eij');
        end
    end
    M_U=sparse(M_U);
end

% Solve the following SDP
%                                min  c'*(x0-barx)
%                                s.t. svec(barQ)+Rq+M_c*c-M_T*vec(T)+M_U*U=sVec(S)
%                                     T>=0, S PSD.
    
% Construct model
prob.c = [x0-barx;zeros((m+1)*m/2+(n+1)*meq,1)]'; % Parameters in the objective function
prob.a = sparse([], [], [], 0, n+(m+1)*m/2+(n+1)*meq); % 0 constraints, n+m*(m+1)/2 scalar variables
prob.blc = []; % Lower bounds for affine constraints
prob.buc = []; % Upper bounds for affine constraints
prob.blx = [-Inf(1,n),zeros(1, (m+1)*m/2),-Inf(1,(n+1)*meq)]; % Lower bounds for scalar variables
prob.bux = Inf(1, n+(m+1)*m/2+(n+1)*meq); % Upper bounds for scalar variables
prob.f = sparse([M_c -M_T M_U;(x0-barx)' zeros(1,(m+1)*m/2+(n+1)*meq)]); % Parameters for scalar variables
prob.g = ([sVec(barQ)+Rq;0])'; % Constant terms in the constraints
prob.accs = [symbcon.MSK_DOMAIN_SVEC_PSD_CONE (n+2)*(n+1)/2 symbcon.MSK_DOMAIN_RPLUS 1]; % PSD cone, Rplus cone


prob.bardim = []; % Dimensions of PSD variables

% Parameters for matrix variables in the objective function
prob.barc.subj = [];
prob.barc.subk = [];
prob.barc.subl = [];
prob.barc.val = [];

% Parameters for matrix variables in the constraints
prob.barf.subi = [];
prob.barf.subj = [];
prob.barf.subk = [];
prob.barf.subl = [];
prob.barf.val = [];


% Parameters setting
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = tol_mosek;
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = tol_mosek;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = tol_mosek;
param.MSK_DPAR_INTPNT_CO_TOL_INFEAS = tol_mosek;
param.MSK_IPAR_AUTO_UPDATE_SOL_INFO = 'MSK_ON';



%[~, res] = mosekopt('minimize echo(5)', prob, param);
[~, res] = mosekopt('minimize echo(0)', prob, param);
sdp_status=0;
val=[];
c=[];
S=[];
if isempty(strfind(res.rcodestr, 'MSK_RES_ERR'))

    solsta = strcat('MSK_SOL_STA_', res.sol.itr.solsta);

    if strcmp(solsta, 'MSK_SOL_STA_OPTIMAL')
        sdp_status=1;
        sol=res.sol.itr.xx;
        val=res.sol.itr.pobjval;
        c=sol(1:n);
        T_sol=sol(n+1:n+(m+1)*m/2);
        T_sol2=T_sol.*(T_sol>=0);
        sol2=[c;T_sol2;sol((m+1)*m/2+n+1:end)];
        %diff=norm(T_sol-T_sol2);
        %fprintf('T difference norm=%4.10f\n\n',diff);
        S=sMat(sVec(barQ)+Rq+[M_c -M_T M_U]*sol2,n+1);

    elseif strcmp(solsta, 'MSK_SOL_STA_DUAL_INFEASIBLE_CER')
        fprintf('Dual infeasibility certificate found.');
        sdp_status=0;
     

    elseif strcmp(solsta, 'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')
        fprintf('Primal infeasibility certificate found.');
        sdp_status=-1;
     
    elseif strcmp(solsta, 'MSK_SOL_STA_UNKNOWN')
        % The solutions status is unknown. The termination code
        % indicates why the optimizer terminated prematurely.
        fprintf('The solution status is unknown.\n');
        fprintf('Termination code: %s (%d) %s.\n', res.rcodestr, res.rcode, res.rmsg);
        sdp_status=-2;
    else
        fprintf('An unexpected solution status is obtained.\n');
    end

else
    fprintf('Error during optimization: %s (%d) %s.\n', res.rcodestr, res.rcode, res.rmsg);
end