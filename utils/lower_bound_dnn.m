function [lb,sdp_status,S,res] = lower_bound_dnn(Q,d,A,b,Aeq,beq,tol_mosek,m,n)

% ======================================================================= %
% Compute the dnn relaxation of the following problem
%     min   x'Qx+2d'x
%     s.t.  Ax <= b
%           Aeq = beq
% INPUT
% Q, d             Parameters of the quadratic objective function
% A, b             Parameters of the inequality constraint A x <= b
% Aeq, beq         Equality constraints: Aeq x = beq
% tol_mosek        Accuracy tolerence used in MOSEK
% m                Size(A,1)
% n                Size(A,2)
%
% OUTPUT
% lb               DNN lower bound given by
%                  lb = max  lambda
%                       s.t. [Q d; d' -lambda]=S+[-A b; zeros(1,n) 1]'T[-A b; zeros(1,n) 1]
%                                + U*[Aeq beq]+[Aeq beq]'*U'
%                            T>=0, S PSD.
% sdp_status
% S                S
% res              res returned by MOSEK
% ======================================================================= %   



[~, res] = mosekopt('symbcon echo(0)');
symbcon = res.symbcon;

barQ=[Q d; d' 0];
barA=[-A b;zeros(1,n) 1];  % barA of dimension (m+1,n+1) 
%                               lb = max  lambda
%                                    s.t.  barQ-[zeros(n,n) zeros(n,1);zeros(1,n) 1]*lambda-barA'*T*barA-U*[Aeq beq]-[Aeq beq]'*U'=S
%                                          T>=0, S PSD.

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


% sVec([zeros(n,n) zeros(n,1);zeros(1,n) 1]*lambda)=M_L*lambda
M_L=zeros(n+1,n+1);
M_L(end,end)=1;
M_L=sVec(M_L);


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

    
% Construct model
prob.c = [-1;zeros((m+1)*m/2+(n+1)*meq,1)]'; % Parameters in the objective function
prob.a = sparse([], [], [], 0, 1+(m+1)*m/2+(n+1)*meq); % 0 constraints, 1+(m+1)*m/2 scalar variables
prob.blc = []; % Lower bounds for affine constraints
prob.buc = []; % Upper bounds for affine constraints
prob.blx = [-Inf,zeros(1, (m+1)*m/2),-Inf(1,(n+1)*meq)]; % Lower bounds for scalar variables
prob.bux = Inf(1, 1+(m+1)*m/2+(n+1)*meq); % Upper bounds for scalar variables
prob.f = sparse([-M_L -M_T M_U]); % Parameters for scalar variables
prob.g = (sVec(barQ))'; % Constant terms in the constraints

prob.accs = [symbcon.MSK_DOMAIN_SVEC_PSD_CONE (n+2)*(n+1)/2]; % PSD cone
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
lb=-Inf;
S=[];
if isempty(strfind(res.rcodestr, 'MSK_RES_ERR'))

    solsta = strcat('MSK_SOL_STA_', res.sol.itr.solsta);

    if strcmp(solsta, 'MSK_SOL_STA_OPTIMAL')
        sdp_status=1;
        lb=-res.sol.itr.pobjval;
        x_sol=res.sol.itr.xx;
        T_sol=x_sol(2:(m+1)*m/2+1);
        T_sol=T_sol.*(T_sol>=0);
        x_sol2=[x_sol(1);T_sol;x_sol((m+1)*m/2+2:end)];
        %diff=norm(x_sol-x_sol2);
        %fprintf('T difference norm=%4.10f\n\n',diff);
        S=sMat(sVec(barQ)+[-M_L -M_T M_U]*x_sol2,n+1);
    

    elseif strcmp(solsta, 'MSK_SOL_STA_DUAL_INFEASIBLE_CER')
        fprintf('Dual infeasibility certificate found.');
        sdp_status=0;
        lb=Inf;

    elseif strcmp(solsta, 'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')
        fprintf('Primal infeasibility certificate found.');
        sdp_status=-1;
     
    elseif strcmp(solsta, 'MSK_SOL_STA_UNKNOWN')
        fprintf('The solution status is unknown.\n');
        fprintf('Termination code: %s (%d) %s.\n', res.rcodestr, res.rcode, res.rmsg);
        sdp_status=-2;
    else
        fprintf('An unexpected solution status is obtained.\n');
    end

else
    fprintf('Error during optimization: %s (%d) %s.\n', res.rcodestr, res.rcode, res.rmsg);
end




end