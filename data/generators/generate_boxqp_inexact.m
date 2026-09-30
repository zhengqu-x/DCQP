function instance = generate_boxqp_inexact(n, seed, options)
%GENERATE_BOXQP_INEXACT Qiu--Yildirim (2024), Algorithm 2.
%   instance = generate_boxqp_inexact(n, seed, options)
%   Objective convention: x'*instance.Q*x + instance.c'*x, 0 <= x <= 1.
%   The paper uses 0.5*x'*Qpaper*x+c'*x, so instance.Q = Qpaper/2.
%   RLT is provably inexact; SDP-RLT need not be inexact. The QP optimizer
%   is NOT planted. A primal/dual certificate of the RLT optimum is saved.
%
%   options.parameter_max (default 10): integer coefficient range 0:M.
%   options.partition (optional): n-vector labelled 1=L, 2=B, 3=U; B nonempty.
%   Otherwise sample each label with equal probability, conditioned on B~=empty,
%   as in Section 6.2. This implements the algorithm, not Julia's random stream.
%   Caller RNG state is restored on return/error. No solvers required.
%
%   Source: Qiu and Yildirim, J Glob Optim 90:293-322 (2024), Algorithm 2
%   and Proposition 16. https://doi.org/10.1007/s10898-024-01407-y

if nargin < 3, options = struct(); end
validateattributes(n, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(seed, {'numeric'}, {'scalar','integer','nonnegative','<=',2^32-1});
if ~isstruct(options) || ~isscalar(options)
    error('dcqp:generationOptions', 'options must be a scalar struct.');
end
allowed = {'parameter_max','partition'};
if ~all(ismember(fieldnames(options), allowed))
    error('dcqp:generationOptions', 'Unknown BoxQP generator option.');
end
M = 10;
if isfield(options, 'parameter_max'), M = options.parameter_max; end
validateattributes(M, {'numeric'}, {'scalar','integer','positive','finite'});
old_rng = rng; restore_rng = onCleanup(@() rng(old_rng));
rng(seed, 'twister');
if isfield(options, 'partition')
    labels = options.partition(:);
    validateattributes(labels, {'numeric'}, {'numel',n,'integer','>=',1,'<=',3});
    if ~any(labels == 2), error('dcqp:generationPartition', 'B must be nonempty.'); end
else
    labels = randi(3,n,1);
    while ~any(labels == 2), labels = randi(3,n,1); end
end
L = find(labels == 1); B = find(labels == 2); U = find(labels == 3);
k = B(randi(numel(B)));
r = zeros(n,1); s = zeros(n,1);
r(U) = randi([0 M],numel(U),1);
s(L) = randi([0 M],numel(L),1);
W = symmetric_random(n,M);
W(L,L) = 0; W(L,B) = 0; W(B,L) = 0; W(k,k) = randi(M);
Y = randi([0 M],n,n);
Y(L,B) = 0; Y(L,U) = 0; Y(B,B) = 0; Y(B,U) = 0;
Z = symmetric_random(n,M);
Z(B,U) = 0; Z(U,B) = 0; Z(U,U) = 0; Z(k,k) = randi(M);
Qpaper = W - Y - Y' + Z;
c = -r + s - W*ones(n,1) + Y'*ones(n,1);
xR = zeros(n,1); xR(B) = 0.5; xR(U) = 1;
XR = zeros(n); XR(B,U) = 0.5; XR(U,B) = 0.5; XR(U,U) = 1;
rlt_primal = 0.5*sum(sum(Qpaper.*XR)) + c'*xR;
rlt_dual = -sum(r) - 0.5*sum(W(:));
assert(abs(rlt_primal-rlt_dual) <= 1e-10*max(1,abs(rlt_dual)), ...
    'dcqp:generationCertificate', 'RLT primal/dual certificate mismatch.');
instance = struct('Q',Qpaper/2,'c',c,'A',[-eye(n);eye(n)], ...
    'b',[zeros(n,1);ones(n,1)],'Aeq',zeros(0,n),'beq',zeros(0,1), ...
    'LB',zeros(n,1),'UB',ones(n,1));
instance.metadata = struct('family','boxqp_inexact','seed',seed, ...
    'source','https://doi.org/10.1007/s10898-024-01407-y', ...
    'construction','Algorithm 2; sampling as in Section 6.2', ...
    'known_optimum',NaN,'known_solution',[], ...
    'rlt_inexact_guaranteed',true,'sdp_rlt_inexact_guaranteed',false, ...
    'known_rlt_optimum',rlt_dual,'parameter_max',M, ...
    'partition',labels,'distinguished_index',k, ...
    'objective_convention','x''*Q*x+c''*x');
instance.certificate = struct('r',r,'s',s,'W',W,'Y',Y,'Z',Z, ...
    'rlt_x',xR,'rlt_X',XR,'rlt_primal',rlt_primal,'rlt_dual',rlt_dual);
end

function S = symmetric_random(n,M)
T = triu(randi([0 M],n,n));
S = T + triu(T,1)';
end
