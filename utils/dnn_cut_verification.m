function [A,b,cut_lb,alpha,action,certified_lb,empty_record] = ...
    dnn_cut_verification(Q,d,A,b,Aeq,beq,c,barx,nu,nuR,delta, ...
    tstar,m,n,independent_validation,konno_added,retry_alpha,p)
%DNN_CUT_VERIFICATION Reproduce the stored DNN cut-validation state machine.
% action is 'retained', 'skipped', or 'terminate'.

cut_lb=Inf;
alpha=0;
action='retained';
certified_lb=-Inf;
empty_record=[];
raw_c_norm=norm(c);
c_norm=max(1,raw_c_norm);

if independent_validation && raw_c_norm==0
    [whole_lb,status,S,~]=lower_bound_dnn( ...
        Q,d,A,b,Aeq,beq,p.mosek_tolerance,m,n,p.mosek_quiet);
    if status~=1
        error('qpsolver:SDPValidation', ...
            'Whole-region SDP validation failed.');
    end
    whole_delta=min(eig(S));
    whole_lb=whole_lb+whole_delta+min(0,whole_delta)*tstar^2;
    certified_lb=whole_lb;
    if whole_lb>=nu
        action='terminate';
    else
        action='skipped';
    end
    return
end

A_discarded=[A;c'/c_norm];
b_discarded=[b;(1+c'*barx)/c_norm];
[~,tstar_cut,status]=gurobilp(-ones(n,1),A_discarded,b_discarded, ...
    Aeq,beq,[],[],p.gurobi_lp_method,p.gurobi_lp_tolerance);
if independent_validation && status==-1
    action='skipped';
    if p.verbose
        fprintf('  SDP cut is redundant: discarded region is empty.\n');
    end
    return
elseif status~=1
    error('qpsolver:SDPValidation','Could not bound the SDP cut region.');
end

cut_lb=nuR+delta+tstar_cut^2*min(delta,0);
if cut_lb>=nu
    % This is the full DNN cut, so its implicit alpha is one. Recording
    % alpha=0 here used to make accepted direct cuts disappear from counts.
    alpha=1;
    A=[A;-c'/c_norm];
    b=[b;(-1-c'*barx)/c_norm];
    [A(end,:),b(end),rescale_status]=rescale_constraint_by_slack( ...
        A(end,:),b(end),A,b,Aeq,beq,p.gurobi_lp_method, ...
        p.gurobi_lp_tolerance);
    retained_empty=(rescale_status==-1);
else
    alpha=1;
    if independent_validation && p.verbose
        fprintf('  Validating SDP cut independently on the current retained region.\n');
    end
    [raw_lb,delta_cut]=discarded_region_bound( ...
        Q,d,A,b,Aeq,beq,c,c_norm,barx,alpha,m,n,p);
    cut_lb=raw_lb+delta_cut+tstar_cut^2*min(delta_cut,0);
    if cut_lb<nu
        alpha=retry_alpha;
        [raw_lb,delta_cut]=discarded_region_bound( ...
            Q,d,A,b,Aeq,beq,c,c_norm,barx,alpha,m,n,p);
        cut_lb=raw_lb+delta_cut+tstar_cut^2*min(delta_cut,0);
    end
    if cut_lb<nu
        if ~p.accept_cut_below_threshold
            error('qpsolver:SDPValidation', ...
                ['DNN cut rejected: corrected discarded-region bound ', ...
                 '%.4e is below nu %.4e.'],cut_lb,nu);
        elseif konno_added
            cut_lb=Inf;
            alpha=0;
            action='skipped';
            if p.verbose
                fprintf('  SDP cut skipped: validated bound remains below nu.\n');
            end
            return
        elseif p.verbose
            fprintf(['  Accepting DNN cut below nu by parameter setting: ', ...
                'corrected bound %.4e < nu %.4e.\n'],cut_lb,nu);
        end
    end

    if independent_validation
        retained_A=[A;-c'/c_norm];
        retained_b=[b;(-alpha-c'*barx)/c_norm];
        [~,~,retained_status]=gurobilp(zeros(n,1),retained_A, ...
            retained_b,Aeq,beq,[],[],p.gurobi_lp_method, ...
            p.gurobi_lp_tolerance);
        if retained_status==-1
            certified_lb=cut_lb;
            % No retained-region inequality is appended in this branch.
            % Keep alpha_record consistent with the number of added DNN cuts.
            alpha=0;
            action='terminate';
            if p.verbose
                fprintf('  Validated SDP cut certifies the remaining region.\n');
            end
            return
        elseif retained_status~=1
            error('qpsolver:SDPValidation', ...
                'Could not check the retained SDP cut region.');
        end
    end

    A=[A;-c'/c_norm];
    b=[b;(-alpha-c'*barx)/c_norm];
    [A(end,:),b(end),rescale_status]=rescale_constraint_by_slack( ...
        A(end,:),b(end),A,b,Aeq,beq,p.gurobi_lp_method, ...
        p.gurobi_lp_tolerance);
    retained_empty=(rescale_status==-1);
end

if retained_empty
    A=A(1:end-1,:);
    b=b(1:end-1);
    [whole_lb,status,S,~]=lower_bound_dnn( ...
        Q,d,A,b,Aeq,beq,p.mosek_tolerance,size(A,1),n,p.mosek_quiet);
    if status~=1
        error('qpsolver:SDPValidation', ...
            'Could not validate the whole region after an empty retained cut.');
    end
    whole_delta=min(eig(S));
    whole_lb=whole_lb+whole_delta+min(whole_delta,0)*tstar^2;
    empty_record=struct('exitflag',-1,'lower_bound',whole_lb, ...
        'target',nu,'validated',whole_lb>=nu,'cut',c,'anchor',barx);
    certified_lb=whole_lb;
    cut_lb=Inf;
    alpha=0;
    if whole_lb>=nu
        action='terminate';
        if p.verbose
            fprintf(['  Retained region is empty: whole-region bound ', ...
                '%.4e >= nu %.4e. Terminating.\n'],whole_lb,nu);
        end
    else
        action='skipped';
        if p.verbose
            fprintf(['  Empty-retained cut skipped: whole-region bound ', ...
                '%.4e < nu %.4e.\n'],whole_lb,nu);
        end
    end
end
end

function [raw_lb,delta] = discarded_region_bound( ...
    Q,d,A,b,Aeq,beq,c,c_norm,barx,alpha,m,n,p)
A_discarded=[A;c'/c_norm];
b_discarded=[b;(alpha+c'*barx)/c_norm];
[raw_lb,status,S,~]=lower_bound_dnn( ...
    Q,d,A_discarded,b_discarded,Aeq,beq,p.mosek_tolerance,m+1,n, ...
    p.mosek_quiet);
if status~=1
    if alpha==1
        message='SDP cut bound was not solved successfully.';
    else
        message='Reduced SDP cut bound was not solved successfully.';
    end
    error('qpsolver:SDPValidation',message);
end
delta=min(eig(S));
if p.verbose
    if alpha==1
        fprintf('  cut lower bound: raw lb = %.4e, min_eig(S) = %.2e\n', ...
            raw_lb,delta);
    else
        fprintf(['  reduced alpha = %.2e: raw cut lb = %.4e, ', ...
            'min_eig(S) = %.2e\n'],alpha,raw_lb,delta);
    end
end
end
