function instance = generate_stqp_many_local(n, seed, options)
%GENERATE_STQP_MANY_LOCAL StQPs with a certified exponential local-minimum count.
%   instance = generate_stqp_many_local(n, seed, options)
%   Objective: x'*instance.Q*x + instance.c'*x on the standard simplex.
%   Pass instance.c/2 as DCQP's d argument. No external solvers are needed.
%
%   This implements the block-combination construction of Bomze, Schachinger
%   and Ullrich, Math. Oper. Res. 43(2):651-674 (2018; online 2017), Section
%   2.2, Eq. (11), Theorem 4, with the seven-dimensional Table 1 seed.
%   https://doi.org/10.1287/moor.2017.0877
%   https://optimization-online.org/wp-content/uploads/2016/05/5452.pdf
%
%   Write n=7*k+r, 0<=r<7. The seven-dimensional seed has 14 strict local
%   maxima, of which seven are global; its two payoff levels are 25/6 and
%   128/31. For r>0, the final block I_r-E_r has r strict global maxima.
%   The returned minimization problem consequently has exactly
%       14^k * max(r,1) strict local minima,
%        7^k * max(r,1) global minima.
%   For n>=7 some strict local minima are nonglobal. For n=6 they are all
%   global. Neither RLT nor SDP-RLT inexactness is claimed for this family.
%
%   options.scale_range    = [0.75 1.25] (positive, sorted two-vector)
%   options.offset_range   = [-0.1 0.1]  (sorted two-vector)
%   options.coupling_margin = 1         (positive scalar)
%   options.permute        = true      (logical scalar)
%   Independent block scales and shifts change their relative importance
%   and the optimal block masses without losing the count guarantee.
%   Permutations also vary variable ordering. No unverified perturbations
%   are added. Caller RNG state is preserved, including after an error.
%
%   Certificate: let R sum coordinates within blocks, A=E_m-m*I_m, and
%   B_i=s_i*S_i+o_i*E. The maximization matrix is
%       G = t*R'*A*R + blkdiag(B_i), and Q=-G.
%   At fixed block masses alpha, maximize each B_i separately. If beta_i
%   is any local payoff and d_i=m*t-beta_i, the reduced objective is
%       t-sum(d_i*alpha_i^2).
%   We choose m*t above every entry of every B_i, so all d_i>0. Its unique
%   maximum is attained at alpha_i=(1/d_i)/sum_j(1/d_j), with all masses
%   positive. Theorem 4 gives exactly the product of the block counts.
%   Substituting each block's global payoff gives the known global value.

if nargin < 3, options = struct(); end
validateattributes(n,{'numeric'},{'scalar','integer','>=',6,'finite'});
validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative','<=',2^32-1,'finite'});
if ~isstruct(options) || ~isscalar(options)
    error('dcqp:generationOptions','options must be a scalar struct.');
end
defaults = struct('scale_range',[0.75 1.25],'offset_range',[-0.1 0.1], ...
    'coupling_margin',1,'permute',true);
unknown = setdiff(fieldnames(options),fieldnames(defaults));
if ~isempty(unknown)
    error('dcqp:generationOptions','Unknown many-local generator option: %s',unknown{1});
end
keys = fieldnames(defaults);
for j = 1:numel(keys)
    if ~isfield(options,keys{j}), options.(keys{j}) = defaults.(keys{j}); end
end
validateattributes(options.scale_range,{'numeric'},{'vector','numel',2,'positive','finite','real'});
validateattributes(options.offset_range,{'numeric'},{'vector','numel',2,'finite','real'});
if options.scale_range(1)>options.scale_range(2) || ...
        options.offset_range(1)>options.offset_range(2)
    error('dcqp:generationOptions','Parameter ranges must be sorted in increasing order.');
end
validateattributes(options.coupling_margin,{'numeric'},{'scalar','positive','finite','real'});
validateattributes(options.permute,{'logical'},{'scalar'});
old_rng = rng; restore_rng = onCleanup(@() rng(old_rng));
rng(seed,'twister');

k = floor(n/7); r = mod(n,7); m = k+(r>0);
block_sizes = [7*ones(1,k), r*ones(1,r>0)];
scales = options.scale_range(1)+diff(options.scale_range(:))*rand(m,1);
offsets = options.offset_range(1)+diff(options.offset_range(:))*rand(m,1);
% Table 1 uses the last-column convention C([1 5 8 8 5 1 0]').
% The equivalent distance convention below puts zeros on the diagonal.
entries = [0 1 5 8 8 5 1];
S = entries(mod((0:6)'-(0:6),7)+1);
p_global = [5;0;2;0;5;0;0]/12;
p_other = [8;8;0;0;15;0;0]/31;
global_payoffs = offsets;
entry_upper_bounds = offsets;
global_payoffs(1:k) = offsets(1:k)+scales(1:k)*(25/6);
entry_upper_bounds(1:k) = offsets(1:k)+8*scales(1:k);
% Use a nonnegative t even when the offsets are negative.
t = (max(0,max(entry_upper_bounds))+options.coupling_margin)/m;
curvature = m*t-global_payoffs;
if any(~isfinite(curvature)) || any(curvature<=0)
    error('dcqp:generationNumerics','Parameters are too large or margin too small for double precision.');
end
weights = (1./curvature)/sum(1./curvature);
G = t*ones(n);
x_global = zeros(n,1); base_block_id = zeros(n,1);
first = 1;
for j = 1:m
    last = first+block_sizes(j)-1; idx = first:last;
    if j<=k
        B = scales(j)*S+offsets(j)*ones(7);
        p = p_global;
    else
        B = scales(j)*(eye(r)-ones(r))+offsets(j)*ones(r);
        p = zeros(r,1); p(1) = 1;
    end
    G(idx,idx) = G(idx,idx)-m*t*ones(block_sizes(j))+B;
    x_global(idx) = weights(j)*p;
    base_block_id(idx) = j;
    first = last+1;
end
Q = -G;
known_optimum = -t+1/sum(1./curvature);
% A concrete strict nonglobal solution is useful for auditing the generator.
x_other = []; other_value = NaN;
if k>0
    other_payoffs = global_payoffs;
    other_payoffs(1) = offsets(1)+scales(1)*(128/31);
    other_curvature = m*t-other_payoffs;
    other_weights = (1./other_curvature)/sum(1./other_curvature);
    x_other = x_global;
    for j = 1:m
        idx = base_block_id==j;
        x_other(idx) = x_global(idx)*(other_weights(j)/weights(j));
    end
    x_other(1:7) = other_weights(1)*p_other;
    other_value = -t+1/sum(1./other_curvature);
end
permutation = 1:n;
if options.permute, permutation = randperm(n); end
Q = Q(permutation,permutation);
x_global = x_global(permutation);
if ~isempty(x_other), x_other = x_other(permutation); end
base_block_id = base_block_id(permutation);
if any(~isfinite(Q(:))) || ~isfinite(known_optimum)
    error('dcqp:generationNumerics','Generated coefficients exceed double precision range.');
end
tol = 1e-10*max(1,norm(Q,'inf'));
assert(abs(x_global'*Q*x_global-known_optimum)<=tol, ...
    'dcqp:generationCertificate','Analytic objective does not match the known solution.');

factor = max(r,1);
instance = struct('Q',Q,'c',zeros(n,1),'A',[-eye(n);eye(n)], ...
    'b',[zeros(n,1);ones(n,1)],'Aeq',ones(1,n),'beq',1, ...
    'LB',zeros(n,1),'UB',ones(n,1));
instance.metadata = struct('family','stqp_many_local','seed',seed, ...
    'source','https://doi.org/10.1287/moor.2017.0877', ...
    'preprint','https://optimization-online.org/wp-content/uploads/2016/05/5452.pdf', ...
    'construction','Section 2.2 Eq. (11), Theorem 4; n=7 seed in Table 1', ...
    'known_optimum',known_optimum,'known_solution',x_global, ...
    'rlt_inexact_guaranteed',false,'sdp_rlt_inexact_guaranteed',false, ...
    'strict_local_minima_count',factor*14^k, ...
    'strict_local_minima_count_exact',integer_power_string(14,k,factor), ...
    'global_minima_count',factor*7^k, ...
    'global_minima_count_exact',integer_power_string(7,k,factor), ...
    'local_count_guaranteed',true,'all_local_minima_global',k==0, ...
    'known_nonglobal_solution',x_other,'known_nonglobal_value',other_value, ...
    'block_sizes',block_sizes,'block_scales',scales,'block_offsets',offsets, ...
    'coupling_parameter',t,'options',options, ...
    'variation','Independent block scales and shifts, then variable permutation', ...
    'objective_convention','x''*Q*x+c''*x');
instance.certificate = struct('base_seven_matrix',S, ...
    'base_global_solution',p_global,'base_global_value',25/6, ...
    'base_nonglobal_solution',p_other,'base_nonglobal_value',128/31, ...
    'block_id',base_block_id,'permutation',permutation, ...
    'block_global_payoffs',global_payoffs,'block_entry_upper_bounds',entry_upper_bounds, ...
    'positive_block_curvature',curvature,'optimal_block_weights',weights);
end

function value = integer_power_string(base, exponent, factor)
% Exact decimal counts without the Symbolic Math Toolbox (numeric fields
% above are rounded doubles once the counts exceed flintmax).
digits = factor;
for j = 1:exponent
    carry = 0;
    for q = 1:numel(digits)
        next = base*digits(q)+carry;
        digits(q) = mod(next,10); carry = floor(next/10);
    end
    while carry>0
        digits(end+1) = mod(carry,10); carry = floor(carry/10); %#ok<AGROW>
    end
end
value = char(fliplr(digits)+'0');
end
