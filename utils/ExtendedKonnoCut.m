function [c,info] = ExtendedKonnoCut(Q,d,A,b,barx,nuR,options)
%EXTENDEDKONNOCUT Construct Proposition 1's cut in the original coordinates.
%   [c,info] = ExtendedKonnoCut(Q,d,A,b,barx,nuR,options)
%
%   Objective: Phi(x) = x'*Q*x + 2*d'*x. All constraints, including bounds,
%   must be supplied as A*x <= b. Eliminate equalities before this call.
%   Q, d and nuR must use the same objective scale and constant convention.
%
%   The retained halfspace is c'*(x-barx) >= 1. The discarded side has
%   Phi(x) >= nuR under Proposition 1's assumptions. This function always
%   sets beta = Phi(barx)-nuR, the LARGEST admissible value (not its negative).
%   It does not modify the solver's constraints or any global settings.
%
%   A nonnegative multiplier representation is found using a simplex LP.
%   Its independent support is extended to a basis G of the active rows.
%   One saddle-point solve computes U from G*U=-I, Q*U+G'*K=0. For each
%   column u_i, an LP minimizes theta_i subject to
%       Q*u_i + theta_i*g + A'*rho_i = 0,
%       (b-A*barx)'*rho_i - beta*theta_i <= lambda_i,
%       theta_i >= 0, rho_i >= 0,
%   where g=Q*barx+d=-G'*lambda. Then c=-G'*theta. These theta_i equal
%   1/mu_i from Proposition 1, with zero for an infeasible mu_i LP.
%   No full constraint basis, inverse, or transformed Hessian is formed.
%
%   Requires Gurobi. OPTIONS fields (unknown fields are rejected):
%     active_tolerance       1e-8  (distance after row normalization)
%     kkt_tolerance          1e-8  (scaled stationarity residual)
%     rank_tolerance         1e-10 (absolute, on normalized active rows)
%     curvature_tolerance    1e-10 (relative to max(1,norm(Q,2)))
%     certificate_tolerance  1e-8
%     lp_tolerance           1e-9  (Gurobi feasibility and optimality)
%     lp_method              1     (0=primal or 1=dual simplex)
%     verbose                false (print construction diagnostics)
%     check_kkt              true  (false: construct a candidate directly;
%                                   independent SDP validation is REQUIRED)
%   With check_kkt=false, use a rank-revealing active basis and solve
%   G'*lambda=-g in least squares to obtain the construction coefficients.
%   Do not test KKT feasibility/stationarity, face curvature, or PSD of the
%   assembled certificate. The input point is unchanged; success then means
%   only that a candidate was computed, not that it is a valid cut.
%
%   On a failed assumption/solve/certificate check: c=[], info.success=false,
%   and info.status/message explain why. Invalid arguments raise an error.
%   On success, info contains beta, theta, mu, basis indices, G, lambda, U,
%   rho and the numerical certificate S,T with W=[-A,b;zeros(1,n),1]:
%     [Q,d;d',-nuR] = S + W'*T*W + (q*h'+h*q')/2,
%     q=[g;beta-g'*barx], h=[-c;1+c'*barx], T>=0, S approximately PSD.
%   Residuals and min(eig(S)) are reported, not hidden. A rigorous bound
%   in finite precision can use the same eigenvalue correction as DCQP.
%   If info.global_bound is true, c=0: the whole feasible set is certified
%   at the target to the reported accuracy. Do not normalize a zero cut.

if nargin<7, options=struct(); end
opts = parse_options(options);
validateattributes(Q,{'double'},{'real','finite','2d','square','nonempty'},mfilename,'Q');
n = size(Q,1);
validateattributes(d,{'double'},{'real','finite','vector','numel',n},mfilename,'d');
validateattributes(barx,{'double'},{'real','finite','vector','numel',n},mfilename,'barx');
validateattributes(A,{'double'},{'real','finite','2d','ncols',n},mfilename,'A');
validateattributes(b,{'double'},{'real','finite','2d'},mfilename,'b');
assert(isempty(b) || isvector(b),'ExtendedKonnoCut:InvalidInput','b must be a vector.');
assert(numel(b)==size(A,1),'ExtendedKonnoCut:InvalidInput','A and b have inconsistent sizes.');
validateattributes(nuR,{'double'},{'real','finite','scalar'},mfilename,'nuR');
Q=full(Q); d=full(d(:)); barx=full(barx(:)); b=full(b(:));
qscale=max(1,norm(Q,2));
assert(norm(Q-Q','fro')<=1e-12*max(1,norm(Q,'fro')), ...
    'ExtendedKonnoCut:NotSymmetric','Q must be symmetric up to roundoff.');
Q=(Q+Q')/2;
vbar=barx'*Q*barx+2*d'*barx;
beta=vbar-nuR;
assert(beta>=0,'ExtendedKonnoCut:InvalidTarget', ...
    'nuR must not exceed Phi(barx). The maximal beta is Phi(barx)-nuR.');
c=[];
info=struct('success',false,'status','initializing','message','', ...
    'vbar',vbar,'nuR',nuR,'beta',beta,'global_bound',false, ...
    'options',opts,'lp_count',0,'coefficient_lp_count',0);
info.kkt_checked=opts.check_kkt;
if opts.verbose
    fprintf('  ExtendedKonnoCut: vbar=%.4e, nu_R=%.4e, beta=%.2e\n',vbar,nuR,beta);
end

% Normalize constraint rows, not variables; this also makes activity tests
% invariant under positive rescaling of the same inequality.
m=size(A,1);
row_scale=full(sqrt(sum(A.^2,2)));
row_scale(row_scale==0)=1;
An=spdiags(1./row_scale,0,m,m)*A;
slack=b-A*barx;
slackn=slack./row_scale;
info.primal_violation=max([0;-slackn]);
if opts.check_kkt && info.primal_violation>opts.active_tolerance
    info=failed(info,'infeasible_point','barx violates the supplied constraints.'); return
end
active=find(abs(slackn)<=opts.active_tolerance & full(sum(An.^2,2))>0);
AI=full(An(active,:));
info.active_indices=active;
g=Q*barx+d;
gradient_scale=max([1,norm(Q*barx,Inf),norm(d,Inf)]);
params=struct('OutputFlag',0,'Method',opts.lp_method, ...
    'FeasibilityTol',opts.lp_tolerance,'OptimalityTol',opts.lp_tolerance);

% Simplex supplies a basic nonnegative representation. Preserve its positive
% support when extending to a full basis: arbitrary QR row selection alone
% need not preserve nonnegative KKT multipliers at a degenerate point.
if opts.check_kkt
lambda_all=zeros(numel(active),1);
if norm(g,Inf)>opts.kkt_tolerance*gradient_scale
    if isempty(active)
        info=failed(info,'not_kkt','There are no active rows to represent the gradient.'); return
    end
    gradient_norm=norm(g,Inf);
    model=struct('A',sparse(AI'),'obj',ones(numel(active),1), ...
        'rhs',full(-g/gradient_norm),'sense',repmat('=',n,1), ...
        'lb',zeros(numel(active),1),'modelsense','min');
    result=solve_lp(model,params);
    info.lp_count=info.lp_count+1;
    info.multiplier_lp_status=result.status;
    if ~strcmp(result.status,'OPTIMAL') || ~isfield(result,'x')
        info=failed(info,'not_kkt',['Nonnegative KKT multiplier LP: ' result.status]); return
    end
    lambda_all=max(0,full(result.x))*gradient_norm;
else
    info.multiplier_lp_status='zero_gradient_within_tolerance';
end
info.kkt_residual=norm(g+AI'*lambda_all,Inf);
if info.kkt_residual>opts.kkt_tolerance*gradient_scale
    info=failed(info,'not_kkt','The KKT stationarity residual exceeds tolerance.'); return
end
support=find(lambda_all>0);
if rank(AI(support,:),opts.rank_tolerance)~=numel(support)
    info=failed(info,'dependent_multiplier_support', ...
        'The simplex multiplier support is numerically dependent.'); return
end
k=rank(AI,opts.rank_tolerance);
basis=support(:);
if ~isempty(AI)
    [~,~,order]=qr(AI','vector');
    for j=order(:)'
        if numel(basis)==k, break; end
        if rank(AI([basis;j],:),opts.rank_tolerance)>numel(basis)
            basis(end+1,1)=j; %#ok<AGROW>
        end
    end
end
G=AI(basis,:);
lambda=lambda_all(basis);
else
    % Construct directly at the supplied point; do not gate on KKT tests.
    k=rank(AI,opts.rank_tolerance);
    if k==0
        basis=zeros(0,1); G=zeros(0,n); lambda=zeros(0,1);
    else
        [~,~,order]=qr(AI','vector');
        basis=order(1:k); basis=basis(:); G=AI(basis,:);
        lambda=-(G'\g);
    end
    info.multiplier_lp_status='not_run';
    info.multiplier_method='active_basis_least_squares';
end
info.k=k;
info.basis_indices=active(basis);
info.basis_row_scale=row_scale(info.basis_indices);
info.G=G;
info.lambda=lambda;
if opts.check_kkt
    info.kkt_residual=norm(g+G'*lambda,Inf);
    if numel(basis)~=k || info.kkt_residual>opts.kkt_tolerance*gradient_scale
        info=failed(info,'basis_failure','Could not construct a suitable active basis.'); return
    end
    Z=null(G,opts.rank_tolerance);
    if isempty(Z)
        info.face_min_eigenvalue=Inf;
    else
        faceQ=Z'*Q*Z; faceQ=(faceQ+faceQ')/2;
        info.face_min_eigenvalue=min(eig(faceQ));
        if info.face_min_eigenvalue<=opts.curvature_tolerance*qscale
            info=failed(info,'face_not_positive_definite', ...
                'Q must be positive definite on the null space of the active rows.'); return
        end
    end
end
if k==0
    U=zeros(n,0);
    info.saddle_residual=0;
else
    saddle=[Q G';G zeros(k)];
    rhs=[zeros(n,k);-eye(k)];
    info.saddle_rcond=rcond(saddle);
    if info.saddle_rcond<eps
        info=failed(info,'singular_saddle_system','The saddle-point system is numerically singular.'); return
    end
    directions=saddle\rhs; % One factorization for all k directions.
    U=directions(1:n,:);
    info.saddle_residual=norm(saddle*directions-rhs,Inf);
    if any(~isfinite(directions(:))) || info.saddle_residual>opts.certificate_tolerance
        info=failed(info,'inaccurate_saddle_solve','The direction equations failed their residual check.'); return
    end
end
info.U=U;
theta=zeros(k,1); rho=zeros(m,k); offset=zeros(k,1);
info.coefficient_lp_status=cell(k,1);
info.lp_equality_residual=zeros(k,1);
info.lp_inequality_residual=zeros(k,1);
for i=1:k
    % The last row is the centered version of the constant-term inequality.
    model=struct('A',sparse([g An';-beta slackn']), ...
        'obj',[1;zeros(m,1)],'rhs',full([-Q*U(:,i);lambda(i)]), ...
        'sense',[repmat('=',n,1);'<'],'lb',zeros(m+1,1),'modelsense','min');
    result=solve_lp(model,params);
    info.lp_count=info.lp_count+1;
    info.coefficient_lp_count=info.coefficient_lp_count+1;
    info.coefficient_lp_status{i}=result.status;
    if ~strcmp(result.status,'OPTIMAL') || ~isfield(result,'x')
        info=failed(info,'coefficient_lp_failed',sprintf('Coefficient LP %d: %s',i,result.status)); return
    end
    values=max(0,full(result.x));
    theta(i)=values(1);
    rho(:,i)=values(2:end)./row_scale;
    offset(i)=lambda(i)+beta*theta(i)-slack'*rho(:,i);
    info.lp_equality_residual(i)=norm(Q*U(:,i)+theta(i)*g+A'*rho(:,i),Inf);
    info.lp_inequality_residual(i)=max(0,-offset(i));
    cert_scale=max([1,norm(Q*U(:,i),Inf),abs(lambda(i))]);
    if any(~isfinite(values)) || ...
            max(info.lp_equality_residual(i),info.lp_inequality_residual(i))>opts.certificate_tolerance*cert_scale
        info=failed(info,'inaccurate_lp_certificate',sprintf('Coefficient LP %d failed its residual check.',i)); return
    end
end
candidate=-G'*theta;
info.theta=theta;
info.mu=Inf(k,1);
positive=theta>0;
info.mu(positive)=1./theta(positive);
info.rho=rho;
info.offset=offset;

% Assemble the certificate in the supplied coordinates. Clipping a tiny
% negative LP offset is accounted for by recomputing S from the identity.
P=sparse(info.basis_indices,(1:k)',1./info.basis_row_scale,m,k);
T=zeros(m+1);
T(1:m,1:m)=full((P*rho'+rho*P')/2);
T(1:m,end)=P*max(offset,0)/2;
T(end,1:m)=T(1:m,end)';
T(1:m+2:end)=0; % Nonnegative square terms are absorbed into S.
W=[-A b;zeros(1,n) 1];
q=[g;beta-g'*barx]; h=[-candidate;1+candidate'*barx];
objective=[Q d;d' -nuR];
S=full(objective-W'*T*W-(q*h'+h*q')/2);
S=(S+S')/2;
info.S=S; info.T=T;
info.certificate_min_eigenvalue=min(eig(S));
info.certificate_identity_residual=norm(objective-S-W'*T*W-(q*h'+h*q')/2,'fro');
info.candidate_c=candidate;
certificate_scale=max(1,norm(objective,2));
if any(~isfinite(S(:))) || (opts.check_kkt && ...
        info.certificate_min_eigenvalue < -opts.certificate_tolerance*certificate_scale)
    info=failed(info,'indefinite_certificate', ...
        'The assembled certificate is not PSD to the requested tolerance.'); return
end
c=candidate;
info.success=true;
info.global_bound=opts.check_kkt && all(theta==0);
if ~opts.check_kkt
    info.status='candidate';
    info.message='Candidate computed without KKT tests; independent SDP validation is required.';
elseif info.global_bound
    info.status='global_bound';
    info.message='Zero cut: the whole feasible set meets the target to the reported certificate accuracy.';
else
    info.status='success';
    info.message='Retain c''*(x-barx)>=1; inspect S for the numerical lower-bound correction.';
end
info.cut_norm=norm(c);
if info.cut_norm>0, info.exclusion_distance=1/info.cut_norm; else, info.exclusion_distance=Inf; end
if opts.verbose
    fprintf('ExtendedKonnoCut: %s, k=%d, ||c||=%.2e, min_eig(S)=%.2e\n', ...
        info.status,k,info.cut_norm,info.certificate_min_eigenvalue);
end
end

function opts=parse_options(options)
assert(isstruct(options) && isscalar(options),'ExtendedKonnoCut:InvalidOptions','options must be a scalar struct.');
opts=struct('active_tolerance',1e-8,'kkt_tolerance',1e-8, ...
    'rank_tolerance',1e-10,'curvature_tolerance',1e-10, ...
    'certificate_tolerance',1e-8,'lp_tolerance',1e-9,'lp_method',1,'verbose',false, ...
    'check_kkt',true);
names=fieldnames(options);
for j=1:numel(names)
    assert(isfield(opts,names{j}),'ExtendedKonnoCut:InvalidOptions','Unknown option: %s',names{j});
    opts.(names{j})=options.(names{j});
end
names={'active_tolerance','kkt_tolerance','rank_tolerance','curvature_tolerance','certificate_tolerance','lp_tolerance'};
for j=1:numel(names)
    validateattributes(opts.(names{j}),{'double'},{'real','finite','scalar','positive'},mfilename,names{j});
end
assert(opts.lp_tolerance>=1e-9 && opts.lp_tolerance<=1e-2, ...
    'ExtendedKonnoCut:InvalidOptions','lp_tolerance must lie in [1e-9,1e-2].');
assert(isscalar(opts.lp_method) && ismember(opts.lp_method,[0 1]), ...
    'ExtendedKonnoCut:InvalidOptions','lp_method must be 0 or 1 (simplex).');
assert(isscalar(opts.verbose) && (islogical(opts.verbose) || ismember(opts.verbose,[0 1])), ...
    'ExtendedKonnoCut:InvalidOptions','verbose must be logical.');
validateattributes(opts.check_kkt,{'logical'},{'scalar'},mfilename,'check_kkt');
end

function result=solve_lp(model,params)
try
    result=gurobi(model,params);
    if strcmp(result.status,'INF_OR_UNBD')
        params.DualReductions=0;
        result=gurobi(model,params);
    end
    if ismember(result.status,{'NUMERIC','SUBOPTIMAL'})
        params.Method=1-params.Method;
        result=gurobi(model,params);
    end
catch exception
    result=struct('status',['ERROR: ' exception.message]);
end
end

function info=failed(info,status,message)
info.status=status; info.message=message;
if info.options.verbose, fprintf('ExtendedKonnoCut: %s: %s\n',status,message); end
end
