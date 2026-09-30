function instance = generate_stqp_inexact(n, seed, options)
%GENERATE_STQP_INEXACT StQP with a planted unique optimum and inexact SDP-RLT.
%   instance = generate_stqp_inexact(n, seed, options)
%   Objective convention: x'*instance.Q*x + instance.c'*x on the simplex.
%   Use d=instance.c/2 when calling DCQP. No solvers are needed to generate
%   the instance or its strict-gap certificate. Caller RNG is restored.
%
%   options.support_size (default min(floor(n/2),n-5)): planted support size,
%       an integer in [2,n-5]; n must be at least 7.
%   options.lambda (default 0): known optimal objective value.
%   options.support_scale (default 0.99): eigenvalue ceiling as a fraction
%       of the separating epsilon; must lie strictly between zero and one.
%
%   Implements the ordinary StQP stage of Algorithm 2, Proposition 5 and
%   Section 4.3(iii) of Bomze, Peng, Qiu and Yildirim (2025), MPC 17:617-651.
%   https://doi.org/10.1007/s12532-025-00285-z
%   The paper subsequently adds a cardinality constraint. This generator
%   intentionally returns the ordinary continuous StQP, before that step.
%   R_AA eigenvalues follow the paper's uniform sampling in (0,0.99*epsilon),
%   rather than matching the different scaling/random stream in Julia code.
%
%   Certificate: R_AB=0, R_AA positive definite, and R_BB copositive imply
%   that xstar uniquely minimizes (x-xstar)'*R*(x-xstar)+lambda on the simplex.
%   The normalized DNN matrix dnn_X is feasible for the StQP SDP-RLT bound
%   and has objective strictly below lambda. Its objective is an UPPER bound
%   on the relaxation optimum, not a lower bound on the original QP.

if nargin < 3, options = struct(); end
validateattributes(n, {'numeric'}, {'real','scalar','integer','>=',7,'finite'});
validateattributes(seed, {'numeric'}, ...
    {'real','scalar','integer','nonnegative','<=',2^32-1,'finite'});
if ~isstruct(options) || ~isscalar(options)
    error('dcqp:generationOptions', 'options must be a scalar struct.');
end
allowed = {'support_size','lambda','support_scale'};
if ~all(ismember(fieldnames(options),allowed))
    error('dcqp:generationOptions', 'Unknown inexact StQP generator option.');
end
k = min(floor(n/2),n-5); lambda = 0; support_scale = 0.99;
if isfield(options,'support_size'), k = options.support_size; end
if isfield(options,'lambda'), lambda = options.lambda; end
if isfield(options,'support_scale'), support_scale = options.support_scale; end
validateattributes(k, {'numeric'}, ...
    {'real','scalar','integer','>=',2,'<=',n-5,'finite'});
validateattributes(lambda, {'numeric'}, {'real','scalar','finite'});
validateattributes(support_scale, {'numeric'}, ...
    {'real','scalar','>',0,'<',1,'finite'});

old_rng = rng; restore_rng = onCleanup(@() rng(old_rng));
rng(seed,'twister');

% Select a random support and normalize positive coordinates. Avoid very
% small coordinates as in the authors' GenerateOptSol.jl implementation.
weights = rand(k,1); weights = weights/sum(weights);
while any(weights < 1e-8)
    weights = rand(k,1); weights = weights/sum(weights);
end
permutation = randperm(n);
support = sort(permutation(1:k));
complement = setdiff(1:n,support);
xstar = zeros(n,1); xstar(support) = weights;

% Equation (53): Horn's exceptional copositive matrix on a 5-cycle.
cycle = diag(ones(4,1),1); cycle(5,1) = 1;
cycle = cycle + cycle';
H = ones(5) - 2*cycle;

% Equation (57), expressed as exact integer ratios before floating-point
% evaluation: diag 7, adjacent entries 4.32, normalized by 78.2.
% F is entrywise nonnegative and positive definite, since its smallest
% eigenvalue is (175-108*(1+sqrt(5))/2)/1955 > 0.
F_numerator = 175*eye(5) + 108*cycle;
F_denominator = 1955;
F = F_numerator/F_denominator;
epsilon = 41/391; % -<H,F>, with sum(F(:))=1.

% Equation (54): an SPD extension and nonnegative off-diagonal block
% preserve copositivity; the Horn principal block preserves exceptionality.
extension_size = n-k-5;
[B, B_eigenvalues] = random_spd(extension_size,3);
C = rand(extension_size,5);
RBB = [B C; C' H];
[RAA, RAA_eigenvalues] = random_spd(k,support_scale*epsilon);
R = zeros(n);
R(support,support) = RAA;
R(complement,complement) = RBB;
P = eye(n) - xstar*ones(1,n);
Q = P'*R*P + lambda*ones(n);
Q = (Q+Q')/2;

% Embed the separating matrix on the five Horn coordinates, as in the
% proof of Proposition 5(ii). X*e supplies the first moment for SDP-RLT.
horn_indices = complement(end-4:end);
dnn_X = zeros(n); dnn_X(horn_indices,horn_indices) = F;
dnn_x = sum(dnn_X,2);
witness_objective = sum(sum(Q.*dnn_X));
support_energy = weights'*RAA*weights;
gap_lower_bound = epsilon-support_energy;
[~, chol_flag] = chol(RAA);
if chol_flag ~= 0 || max(RAA_eigenvalues) >= epsilon || ...
        gap_lower_bound <= 0 || witness_objective >= lambda
    error('dcqp:generationCertificate', ...
        'Strict-gap certificate lost numerically; change the options or seed.');
end

instance = struct('Q',Q,'c',zeros(n,1),'A',[-eye(n);eye(n)], ...
    'b',[zeros(n,1);ones(n,1)],'Aeq',ones(1,n),'beq',1, ...
    'LB',zeros(n,1),'UB',ones(n,1));
instance.metadata = struct('family','stqp_inexact','seed',seed, ...
    'source','https://doi.org/10.1007/s12532-025-00285-z', ...
    'source_code',['https://github.com/newfound21/' ...
        'Tractable_relaxations_for_sparse_StQP/blob/main/src/Generate_Q_COP.jl'], ...
    'construction','Algorithm 2 ordinary StQP stage; Proposition 5; Section 4.3(iii)', ...
    'known_optimum',lambda,'known_solution',xstar,'unique_optimum',true, ...
    'rlt_inexact_guaranteed',true,'sdp_rlt_inexact_guaranteed',true, ...
    'support_size',k,'support_scale',support_scale, ...
    'cardinality_constraint_added',false, ...
    'objective_convention','x''*Q*x+c''*x');
instance.certificate = struct('support',support,'complement',complement, ...
    'R',R,'RAA',RAA,'RBB',RBB,'extension_B',B,'extension_C',C, ...
    'RAA_eigenvalues',RAA_eigenvalues,'B_eigenvalues',B_eigenvalues, ...
    'horn_indices',horn_indices,'horn_matrix',H, ...
    'F_numerator',F_numerator,'F_denominator',F_denominator, ...
    'epsilon',epsilon,'dnn_x',dnn_x,'dnn_X',dnn_X, ...
    'witness_objective',witness_objective, ...
    'witness_objective_analytic',lambda+support_energy-epsilon, ...
    'gap_lower_bound',gap_lower_bound, ...
    'known_solution_objective',xstar'*Q*xstar);
end

function [S, eigenvalues] = random_spd(n,ceiling)
% Orthogonal conjugation of positive, uniformly sampled eigenvalues.
if n == 0, S = zeros(0); eigenvalues = zeros(0,1); return; end
eigenvalues = ceiling*max(rand(n,1),eps);
[V,~] = qr(randn(n));
S = V*diag(eigenvalues)*V';
S = (S+S')/2;
end
