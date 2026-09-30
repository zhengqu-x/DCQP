function manifest = prepare_concave_boxqp_unique(n)
%PREPARE_CONCAVE_BOXQP_UNIQUE Save one strictly concave BoxQP with a unique optimum.
%   manifest = prepare_concave_boxqp_unique(n)
%
%   Uses the odd-dimensional SDP-RLT-inexact example in Qiu--Yildirim
%   (JOGO 2024, Section 6.3), plus two explicit perturbations. An affine
%   tie-breaker makes the planted binary optimizer unique. The term
%   epsilon*sum(x_i*(1-x_i)) makes the objective strictly concave while
%   preserving both the vertex optimum and the SDP-RLT witness value.
%   For even n, the active construction has dimension n-1 and the last
%   variable is pinned to zero by the tie-breaker. The SDP-RLT gap remains
%   provably positive for every integer n >= 3.
%
%   The saved objective convention is x'*Q*x + c'*x, with no constant.
%   Adding metadata.objective_offset recovers the unshifted construction.
%   The returned manifest and instance are saved together under
%   data/structured/generated/concave_boxqp_unique/. Existing, different
%   files are not overwritten.
%
%   Source: https://doi.org/10.1007/s10898-024-01407-y

validateattributes(n,{'numeric'},{'scalar','integer','>=',3,'finite'});
active_n = n - double(mod(n,2)==0);
k = (active_n-1)/2;
epsilon = 1e-3;
delta = 1/(8*active_n^2);

v = zeros(n,1);
v(1:k) = 1;
e = ones(n,1);
Qpaper = zeros(n);
Qpaper(1:active_n,1:active_n) = ones(active_n)/active_n-eye(active_n);
Q = Qpaper/2-epsilon*eye(n);
c = epsilon*e+delta*(e-2*v);
objective_offset = delta*k;

% This is the SDP-RLT witness from Qiu--Yildirim, extended by a zero row
% and column when n is even. Its diagonal equals x, so the epsilon term
% vanishes in the lifted objective, just as it does at every vertex.
x_sdp = zeros(n,1);
x_sdp(1:active_n) = 0.5;
X_sdp = zeros(n);
X_sdp(1:active_n,1:active_n) = ...
    0.25*ones(active_n)+(active_n*eye(active_n)-ones(active_n))/(4*(active_n-1));
known_optimum = v'*Q*v+c'*v;
sdp_witness_value = sum(sum(Q.*X_sdp))+c'*x_sdp;
proven_gap = known_optimum-sdp_witness_value;
expected_gap = 1/(16*active_n);
tol = 1e-10;
assert(abs(known_optimum+objective_offset+(active_n^2-1)/(8*active_n))<tol, ...
    'dcqp:concaveBoxQP','Planted objective mismatch.');
assert(abs(proven_gap-expected_gap)<tol && proven_gap>0, ...
    'dcqp:concaveBoxQP','SDP-RLT gap certificate mismatch.');
assert(max(eig(Q))<0 && min(eig(X_sdp-x_sdp*x_sdp'))>-tol, ...
    'dcqp:concaveBoxQP','Concavity or SDP certificate failed.');
assert(all(diag(X_sdp)<=x_sdp+tol) && min(X_sdp(:))>=-tol && ...
    min(reshape(x_sdp*e'-X_sdp,[],1))>=-tol && ...
    min(reshape(X_sdp-x_sdp*e'-e*x_sdp'+e*e',[],1))>=-tol, ...
    'dcqp:concaveBoxQP','SDP witness violates a box RLT inequality.');

instance = struct('Q',Q,'c',c,'A',[-eye(n);eye(n)], ...
    'b',[zeros(n,1);e],'Aeq',zeros(0,n),'beq',zeros(0,1), ...
    'LB',zeros(n,1),'UB',e);
instance.metadata = struct('family','concave_boxqp_unique', ...
    'source','https://doi.org/10.1007/s10898-024-01407-y', ...
    'construction','Qiu--Yildirim Section 6.3 with strict-concavity and affine tie-breaker perturbations', ...
    'active_n',active_n,'known_solution',v,'known_optimum',known_optimum, ...
    'unshifted_known_optimum',known_optimum+objective_offset, ...
    'objective_offset',objective_offset,'epsilon',epsilon,'delta',delta, ...
    'strictly_concave',true,'unique_global_solution',true, ...
    'sdp_rlt_inexact_guaranteed',true,'proven_sdp_rlt_gap',proven_gap, ...
    'objective_convention','x''*Q*x+c''*x');
instance.certificate = struct('sdp_x',x_sdp,'sdp_X',X_sdp, ...
    'sdp_witness_value',sdp_witness_value, ...
    'known_optimum',known_optimum,'proven_gap',proven_gap);

data_dir = fileparts(fileparts(mfilename('fullpath')));
folder = fullfile(data_dir,'structured','generated','concave_boxqp_unique');
if ~isfolder(folder), mkdir(folder); end
id = sprintf('concave_boxqp_unique_n%d',n);
file = fullfile(folder,[id '.mat']);
manifest = struct('id',id,'family','concave_boxqp_unique', ...
    'n',n,'seed',0,'file',file);
if isfile(file)
    previous = load(file,'instance','manifest');
    if ~isfield(previous,'instance') || ~isfield(previous,'manifest') || ...
            ~isequaln(previous.instance,instance) || ~isequaln(previous.manifest,manifest)
        error('dcqp:generationOverwrite','Different instance already exists at %s.',file);
    end
else
    save(file,'instance','manifest');
end
fprintf('Ready: %s\n',file);
fprintf('Known optimum %.12g; SDP-RLT gap at least %.12g.\n', ...
    known_optimum,proven_gap);
end
