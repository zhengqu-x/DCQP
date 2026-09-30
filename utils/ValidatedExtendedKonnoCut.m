function [c,info] = ValidatedExtendedKonnoCut(Q,d,A,b,barx,nuR,r0,options)
%VALIDATEDEXTENDEDKONNOCUT Validate and enlarge an Extended Konno cut.
%   [c,info] = ValidatedExtendedKonnoCut(Q,d,A,b,barx,nuR,r0,options)
%   First construct c0 with ExtendedKonnoCut, using beta=Phi(barx)-nuR.
%   Independently solve the DNN bound on each discarded region
%       A*x<=b, (c0/10^k)'*(x-barx)<=1, k=0,1,2,...
%   Accept only if raw_lb+delta+min(0,delta)*r0^2 >= nuR.
%   Stop at the first failed check or a validated norm(c)<1. Return the
%   last validated c, or [] if none passed. A failed solve is not a proof
%   that the cut is invalid; it stops this certification procedure.
%
%   Objective: x'*Q*x+2*d'*x. All constraints, including bounds, must be
%   in A,b; eliminate equalities first. The caller must supply a VALID
%   upper bound r0 on norm(x) over A*x<=b (e.g. 1 on the unit simplex).
%   The retained halfspace is c'*(x-barx)>=1. This function does not
%   append any constraints to the DCQP solver.
%
%   OPTIONS: sdp_tolerance=1e-8, konno_options=struct(), verbose=true,
%   mosek_quiet=false.
%   Requires Gurobi, MOSEK, and utils/lower_bound_dnn on the MATLAB path.
%   info.history records every attempted factor, including a failed one;
%   info.certificates stores its reconstructed S and raw solver result.
%   These are numerical certificates, using DCQP's eigenvalue correction.

if nargin<8, options=struct(); end
opts=struct('sdp_tolerance',1e-8,'konno_options',struct(),'verbose',true, ...
    'mosek_quiet',false);
assert(isstruct(options) && isscalar(options), ...
    'ValidatedExtendedKonnoCut:InvalidOptions','options must be a scalar struct.');
names=fieldnames(options);
for j=1:numel(names)
    assert(isfield(opts,names{j}),'ValidatedExtendedKonnoCut:InvalidOptions', ...
        'Unknown option: %s',names{j});
    opts.(names{j})=options.(names{j});
end
validateattributes(r0,{'double'},{'real','finite','scalar','nonnegative'});
validateattributes(opts.sdp_tolerance,{'double'},{'real','finite','scalar','positive'});
validateattributes(opts.verbose,{'logical'},{'scalar'});
validateattributes(opts.mosek_quiet,{'logical'},{'scalar'});
[c0,konno]=ExtendedKonnoCut(Q,d,A,b,barx,nuR,opts.konno_options);
c=[];
info=struct('success',false,'status','construction_failed','message',konno.message, ...
    'konno',konno,'initial_cut',c0,'radius_bound',r0,'target',nuR, ...
    'best_factor',NaN,'best_norm',NaN,'best_lower_bound',NaN, ...
    'global_bound',false,'history',struct([]),'certificates',{{}},'options',opts);
if ~konno.success, return; end
barx=full(barx(:)); d=full(d(:)); b=full(b(:)); Q=full((Q+Q')/2);
assert(norm(barx)<=r0+1e-10*max(1,r0), ...
    'ValidatedExtendedKonnoCut:InvalidRadius','r0 does not even contain barx.');
validateattributes(norm(c0),{'double'},{'finite'});
n=numel(barx); factor=1;
while true
    candidate=c0/factor;
    cut_norm=norm(candidate);
    if cut_norm==0
        A_cut=A; b_cut=b; % Zero cut: certify the entire feasible region.
    else
        A_cut=[A;candidate'/cut_norm];
        b_cut=[b;(1+candidate'*barx)/cut_norm];
    end
    row=struct('factor',factor,'norm_c',cut_norm,'status','not_solved', ...
        'raw_lower_bound',NaN,'min_eigenvalue',NaN,'correction',NaN, ...
        'lower_bound',NaN,'passes',false,'seconds',NaN,'message','');
    S=[]; res=struct(); started=tic;
    try
        [raw,status,S,res]=lower_bound_dnn(Q,d,A_cut,b_cut,zeros(0,n), ...
            zeros(0,1),opts.sdp_tolerance,size(A_cut,1),n,opts.mosek_quiet);
        if isfield(res,'sol') && isfield(res.sol,'itr') && isfield(res.sol.itr,'solsta')
            row.status=res.sol.itr.solsta;
        end
        if status==1 && isfinite(raw) && ~isempty(S) && all(isfinite(S(:)))
            delta=min(eig((S+S')/2));
            row.raw_lower_bound=raw; row.min_eigenvalue=delta;
            row.correction=delta+min(0,delta)*r0^2;
            row.lower_bound=raw+row.correction;
            row.passes=isfinite(row.lower_bound) && row.lower_bound>=nuR;
        end
    catch exception
        row.status='error'; row.message=exception.message;
    end
    row.seconds=toc(started);
    if isempty(info.history), info.history=row; else, info.history(end+1)=row; end
    info.certificates{end+1}=struct('S',S,'solver_result',res);
    if opts.verbose
        fprintf('  Konno validation: factor=%.2e, norm(c)=%.2e, lb=%.4e, target=%.4e, pass=%d\n', ...
            factor,cut_norm,row.lower_bound,nuR,row.passes);
    end
    if ~row.passes
        if isempty(c)
            info.status='initial_validation_failed';
            info.message='The original Konno cut did not pass independent DNN validation.';
        else
            info.status='enlargement_validation_failed';
            info.message='Returning the last validated cut; the next enlargement did not pass.';
        end
        return
    end
    c=candidate; info.success=true;
    info.best_factor=factor; info.best_norm=cut_norm; info.best_lower_bound=row.lower_bound;
    info.global_bound=(cut_norm==0);
    if cut_norm<1
        info.status='norm_below_one';
        info.message='Returning a validated cut with norm strictly below one.';
        return
    end
    if ~isfinite(factor*10)
        info.status='factor_overflow'; info.message='Returning the last validated cut before factor overflow.';
        return
    end
    factor=factor*10;
end
end
