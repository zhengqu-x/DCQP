function [bestub,best_sol,bestlb,nb_iters,cut_counts,diagnostics]=qpsolver(Q,d,A,b,Aeq,beq,lb,ub,sol,parameters)

%==========================================================================%
% Solve the following quadratic program:
% min  x'*Q*x + 2*d'*x
% s.t. Ax <= b
%      Aeq x = beq
%
% INPUT
% Q, d             Parameters of the quadratic objective function
% A, b             Parameters of the inequality constraint A x <= b
% Aeq, beq         Equality constraints: Aeq x = beq
% lb, ub           Known lower and upper bounds of the objective value
% sol              Initial feasible solution
% parameters       Structure containing algorithm parameters:
%
% OUTPUT
% bestub           Best upper bound of the objective value found
% best_sol         Corresponding solution vector achieving bestub
% bestlb           Best lower bound of the objective value found
% nb_iters         Number of outer iterations performed
% cut_counts       Added-cut counts: total, dnn, and konno
% diagnostics      Detailed iteration records for the caller to save
%==========================================================================

met_glp=parameters.gurobi_lp_method;
max_iterations=parameters.max_iterations;
max_time=parameters.max_time;
spn=parameters.dc_regularization;
gap_tol=parameters.gap_tolerance;
tol_mosek=parameters.mosek_tolerance;
tol_glp=parameters.gurobi_lp_tolerance;
eta=parameters.eta;
known_start=isfield(parameters,'known_solution') && ~isempty(parameters.known_solution);
konnofirst=parameters.konnofirst;


bestub= ub;
bestlb= lb;


lb_record=zeros(max_iterations,1);
ub_record=zeros(max_iterations,1);
bestub_record=zeros(max_iterations,1);
bestlb_record=zeros(max_iterations,1);
cut_lb_record=Inf(max_iterations,1);
alpha_record=zeros(max_iterations,1);
konno_lb_record=Inf(max_iterations,1);
konno_record=cell(max_iterations,1);
time_record_konno=zeros(max_iterations,1);
konno_recovery_record=cell(max_iterations,1);
time_record_konno_recovery=zeros(max_iterations,1);
dc_recovery_record=cell(max_iterations,1);
time_record_dc_recovery=zeros(max_iterations,1);
empty_retained_record=cell(max_iterations,1);

time_record_lb=zeros(max_iterations,1);
time_record_cut_lb=zeros(max_iterations,1);
time_record_kkt=zeros(max_iterations,1);
time_record_generate_cut=zeros(max_iterations,1);




A_bar=A;
b_bar=b;
initial_inequality_count=size(A_bar,1);



cut_val=Inf;

n=size(A_bar,2);
[~,tstar0,~]=gurobilp(-ones(n,1),A_bar,b_bar,Aeq,beq,[],[],met_glp,tol_glp);

solution_scale=1;
if abs(tstar0)>1
    solution_scale=ceil(abs(tstar0));
    A_bar=A_bar*solution_scale;
    Aeq=Aeq*solution_scale;
    Q=Q*solution_scale^2;
    d=d*solution_scale;
end
best_sol=sol/solution_scale;

[M,N] = DC_decomposition(Q,spn);


for j=1:size(A_bar,1)
    [A_bar(j,:),b_bar(j),rescale_status] = rescale_constraint_by_slack(A_bar(j,:),b_bar(j),A_bar,b_bar,Aeq,beq,met_glp,tol_glp);
    if rescale_status==-1
        error('qpsolver:EmptyInputRegion','The input region is empty during initial constraint rescaling.');
    end
end



i=0;
if isfield(parameters,'solve_timer')
    solve_timer=parameters.solve_timer;
else
    solve_timer=tic;
end
while i<max_iterations && toc(solve_timer)<max_time
    i=i+1;

    m=size(A_bar,1);
    n=size(A_bar,2);

    tic

    [~,tstar,exitflag]=gurobilp(-ones(n,1),A_bar,b_bar,Aeq,beq,[],[],met_glp,tol_glp);
    if exitflag==-1
        if i==1
            error('The feasible region is empty.\n');
        else
            lb_record(i)=min(cut_lb_record(1:i-1));
            break
        end
    end

    if parameters.verbose
        fprintf('iteration %d\n',i);
    end

    if cut_val<1 || mod(i-1,10)==0
        [lb_record(i),barz]=compute_lower_bound( ...
            Q,d,A_bar,b_bar,Aeq,beq,tol_mosek,m,n,tstar, ...
            parameters.verbose,parameters.mosek_quiet);
    else
        lb_record(i)=lb_record(i-1);
    end


    bestlb=max(bestlb,lb_record(i));
    bestlb_record(i)=bestlb;
    if parameters.verbose==true
        fprintf('  corrected lb = %.4e, bestlb = %.4e\n', ...
            lb_record(i),bestlb);
    end

    time_record_lb(i)=toc;

    tic
    [barx,v]=compute_upper_bound( ...
        Q,d,A_bar,b_bar,Aeq,beq,M,N,barz,best_sol,bestub,i,known_start,parameters);
    rho=linearization_minimum(Q,d,barx,A_bar,b_bar,Aeq,beq,met_glp,tol_glp);

    time_record_kkt(i)=toc;


    if v<bestub
        best_sol=barx;
    end
    bestub=min(v,bestub);
    ub_record(i)=v;
    bestub_record(i)=bestub;
    if parameters.verbose==true
        if abs(bestub)>gap_tol
            displayed_relative_gap=abs(bestub-bestlb)/abs(bestub);
        else
            displayed_relative_gap=abs(bestub-bestlb);
        end
        fprintf(['  v = %.4e, bestub = %.4e, bestlb = %.4e, ', ...
            'relative gap = %.4e\n'], ...
            v,bestub,bestlb,displayed_relative_gap);
    end


    tic;
    if bestub<=bestlb+abs(bestub)*gap_tol
        break
    end

    [nu,nuR,beta]=dcqp_cut_targets(bestub,v,eta,gap_tol);


    konno_added=false;
    konno_refined=false;
    konno_threshold=max(1e-5,abs(bestub)*gap_tol*0.1);
    if konnofirst && abs(bestub-v)<=konno_threshold
        konno_timer=tic;
        [A_bar,b_bar,barx,v,bestub,best_sol,nu,nuR,beta,konno]= ...
            compute_konno_cut(Q,d,A_bar,b_bar,Aeq,beq,M,N,barx,v, ...
            bestub,best_sol,nu,nuR,beta,barz,tstar,parameters);
        konno_added=konno.added;
        konno_refined=konno.refined;
        konno_record{i}=konno.record;
        konno_recovery_record{i}=konno.recovery;
        time_record_konno_recovery(i)=konno.recovery_time;
        if konno_refined
            ub_record(i)=v;
            bestub_record(i)=bestub;
        end
        time_record_konno(i)=toc(konno_timer);
        if konno_record{i}.accepted
            konno_lb_record(i)=konno_record{i}.lower_bound;
            if parameters.verbose
                fprintf('  Konno accepted: factor=%.2e, norm(c)=%.2e, validated lb=%.4e\n', ...
                    konno_record{i}.factor,konno_record{i}.norm_c,konno_lb_record(i));
            end
            if konno_record{i}.region_empty
                bestlb=max(bestlb,konno_lb_record(i));
                lb_record(i)=bestlb; bestlb_record(i)=bestlb;
                if parameters.verbose, fprintf('  Konno certifies the entire remaining region.\n'); end
                break
            end
            m=size(A_bar,1); % Include the newly retained halfspace in the SDP.
            % Keep the same anchor and relaxation point across the two cuts.
            konno_record{i}.sdp_region_rows=m;
            konno_record{i}.sdp_objective_point=barz;
            if parameters.verbose
                fprintf(['  DNN cut after generalized Konno cut: retained region ', ...
                    'has %d inequalities; keeping the same barx and barz.\n'],m);
            end
        elseif parameters.verbose
            fprintf('  Konno skipped: %s; %s\n',konno_record{i}.status,konno_record{i}.message);
        end
    elseif konnofirst && parameters.verbose
        fprintf('  Konno skipped: objective is outside the closeness threshold.\n');
    end
    rho=linearization_minimum(Q,d,barx,A_bar,b_bar,Aeq,beq,met_glp,tol_glp);
    % Roll back only the attempted SDP cut if its generation/validation fails.
    A_before_sdp=A_bar; b_before_sdp=b_bar;
    generation_timer=tic; % Separate Konno time from SDP generation time.
    try

    tol_mosek_cut=parameters.tol_mosek_cut;
    if parameters.verbose==true
        % Print before the solve so these values remain visible if MOSEK fails.
        fprintf('  nu = %.4e, nu_R = %.4e\n',nu,nuR);
        fprintf('  beta = %.2e, v - nu_R = %.4e, nu_R - nu = %.4e\n', ...
            beta, v-nuR, nuR-nu);
        fprintf('  rho = %.2e, min(rho+beta,0) = %.2e\n',rho,min(rho+beta,0));
    end
    [c,cut_val,~,S] = generate_cut_dnn(Q,d,A_bar,b_bar,Aeq,beq,m,n, ...
        nuR,barx,barz,tol_mosek_cut,beta,rho,parameters.mosek_quiet);

    while tol_mosek_cut<=1e-5 && isempty(c)
        tol_mosek_cut=tol_mosek_cut*10;
        fprintf('  retrying cut generation with tol_mosek_cut = %.2e\n',tol_mosek_cut);
        [c,cut_val,~,S] = generate_cut_dnn(Q,d,A_bar,b_bar,Aeq,beq,m,n, ...
            nuR,barx,barz,tol_mosek_cut,beta,rho,parameters.mosek_quiet);
    end
    dc_refined=false;
    if isempty(c) && ~konno_added
        recovery_timer=tic;
        [recovered_point,dc_recovery_record{i}]=refine_failed_cut_point( ...
            Q,d,A_bar,b_bar,Aeq,beq,M,N,barx,beta,parameters);
        time_record_dc_recovery(i)=toc(recovery_timer);
        dc_recovery_record{i}.sdp_objective_point=barz;
        dc_recovery_record{i}.konno_already_added=konno_added;
        if dc_recovery_record{i}.refined
            % This is a new recovery attempt; keep the existing SDP point
            % and all previously validated cuts on the retained region.
            dc_refined=true;
            barx=recovered_point;
            v=barx'*Q*barx+2*d'*barx;
            if v<bestub, best_sol=barx; end
            bestub=min(bestub,v);
            ub_record(i)=v; bestub_record(i)=bestub;
            % Apply the unchanged cut-target/beta policy to the new anchor.
            [nu,nuR,beta]=dcqp_cut_targets(bestub,v,eta,gap_tol);
            rho=linearization_minimum(Q,d,barx,A_bar,b_bar,Aeq,beq,met_glp,tol_glp);
            dc_recovery_record{i}.beta_after=beta;
            dc_recovery_record{i}.nu_after=nu;
            dc_recovery_record{i}.nuR_after=nuR;
            if parameters.verbose
                fprintf('  SDP recovery retry: bestub=%.4e, v=%.4e, nu=%.4e, nu_R=%.4e, beta=%.2e; barz unchanged.\n', ...
                    bestub,v,nu,nuR,beta);
            end
            tol_mosek_cut=parameters.tol_mosek_cut;
            [c,cut_val,~,S]=generate_cut_dnn(Q,d,A_bar,b_bar,Aeq,beq,m,n, ...
                nuR,barx,barz,tol_mosek_cut,beta,rho,parameters.mosek_quiet);
            while tol_mosek_cut<=1e-5 && isempty(c)
                tol_mosek_cut=tol_mosek_cut*10;
                fprintf('  retrying recovered cut with tol_mosek_cut = %.2e\n',tol_mosek_cut);
                [c,cut_val,~,S]=generate_cut_dnn(Q,d,A_bar,b_bar,Aeq,beq,m,n, ...
                    nuR,barx,barz,tol_mosek_cut,beta,rho,parameters.mosek_quiet);
            end
            dc_recovery_record{i}.retry_success=~isempty(c);
        end
    end
    independent_validation=konno_added || konno_refined || dc_refined;
    if isempty(c)
        error('qpsolver:SDPCutGeneration','Failed to generate cut with MOSEK tolerance %4.2e. Try decreasing the error tolerance epsilon.',tol_mosek_cut)
    else
        delta=min(eig(S));
        if parameters.verbose==true
            fprintf('  cut result: min_eig(S) = %.2e, cut_val = %.2e\n', ...
                delta,cut_val);
        end
    end


    time_record_generate_cut(i)=toc(generation_timer);

    verification_timer=tic;
    [A_bar,b_bar,cut_lb_record(i),alpha_record(i),action, ...
        certified_lb,empty_retained_record{i}]=dnn_cut_verification( ...
        Q,d,A_bar,b_bar,Aeq,beq,c,barx,nu,nuR,delta,tstar,m,n, ...
        independent_validation,konno_added,0.9,parameters);

    if isfinite(certified_lb)
        bestlb=max(bestlb,certified_lb);
        lb_record(i)=bestlb;
        bestlb_record(i)=bestlb;
    end
    if strcmp(action,'terminate')
        time_record_cut_lb(i)=toc(verification_timer);
        break
    elseif strcmp(action,'skipped')
        cut_val=0;
        time_record_cut_lb(i)=toc(verification_timer);
        continue
    end

    % A Konno cut changed the region: refresh the main relaxation next time.
    if konno_added, cut_val=0; end


    if parameters.verbose==true
        fprintf('  cut retained: cut_lb = %.4e, delta = %.2e\n', ...
            cut_lb_record(i),delta);
    end
    time_record_cut_lb(i)=toc(verification_timer);
    catch cut_error
        if ~ismember(cut_error.identifier,{'qpsolver:SDPCutGeneration','qpsolver:SDPValidation'})
            rethrow(cut_error)
        end
        A_bar=A_before_sdp; b_bar=b_before_sdp;
        rejected_cut_val=cut_val;
        cut_lb_record(i)=Inf; alpha_record(i)=0; cut_val=0;
        if time_record_generate_cut(i)==0
            time_record_generate_cut(i)=toc(generation_timer);
        end
        if ~konno_added
            % SDP recovery above checked rho and, when needed, refined barx.
            % Validation failures also need that check before forced Konno.
            fallback_timer=tic;
            state=struct('barx',barx,'best_sol',best_sol,'bestub',bestub, ...
                'v',v,'nu',nu,'nuR',nuR,'beta',beta);
            previous_konno=konno_record{i};
            [A_bar,b_bar,state,fallback_record,fallback_recovery]=dcqp_konno_fallback( ...
                Q,d,A_bar,b_bar,Aeq,beq,M,N,state,barz,abs(tstar),parameters, ...
                isempty(dc_recovery_record{i}),isempty(konno_recovery_record{i}));
            time_record_konno(i)=time_record_konno(i)+toc(fallback_timer);
            best_sol=state.best_sol; bestub=state.bestub; v=state.v;
            ub_record(i)=v; bestub_record(i)=bestub;
            fallback_record.previous_attempt=previous_konno;
            fallback_record.sdp_error_identifier=cut_error.identifier;
            fallback_record.sdp_error_message=cut_error.message;
            konno_record{i}=fallback_record;
            if ~isempty(fallback_recovery)
                if isempty(konno_recovery_record{i})
                    konno_recovery_record{i}=fallback_recovery;
                else
                    konno_recovery_record{i}.sdp_fallback_recovery=fallback_recovery;
                end
            end
            if fallback_record.accepted
                konno_lb_record(i)=fallback_record.lower_bound;
                if parameters.verbose
                    fprintf('  Konno accepted: factor=%.2e, norm(c)=%.2e, validated lb=%.4e (SDP failure fallback)\n', ...
                        fallback_record.factor,fallback_record.norm_c,konno_lb_record(i));
                end
                if fallback_record.region_empty
                    bestlb=max(bestlb,konno_lb_record(i));
                    lb_record(i)=bestlb; bestlb_record(i)=bestlb;
                    if parameters.verbose, fprintf('  Fallback Konno certifies the entire remaining region.\n'); end
                    break
                end
                konno_added=true;
            end
        end
        if ~konno_added
            if strcmp(cut_error.identifier,'qpsolver:SDPValidation') && ~isempty(c)
                final_parameters=parameters;
                final_parameters.accept_cut_below_threshold=true;
                if parameters.verbose
                    fprintf(['  Generalized Konno fallback failed; retrying ', ...
                        'the DNN cut with alpha = 0.5.\n']);
                end
                verification_timer=tic;
                [A_bar,b_bar,cut_lb_record(i),alpha_record(i),action, ...
                    certified_lb,empty_retained_record{i}]=dnn_cut_verification( ...
                    Q,d,A_before_sdp,b_before_sdp,Aeq,beq,c,barx,nu,nuR, ...
                    delta,tstar,m,n,independent_validation,false,0.5, ...
                    final_parameters);
                time_record_cut_lb(i)=toc(verification_timer);
                cut_val=rejected_cut_val;
                if isfinite(certified_lb)
                    bestlb=max(bestlb,certified_lb);
                    lb_record(i)=bestlb;
                    bestlb_record(i)=bestlb;
                end
                if strcmp(action,'terminate')
                    break
                elseif strcmp(action,'skipped')
                    cut_val=0;
                    continue
                end
                if parameters.verbose
                    fprintf(['  Accepted final DNN retry: alpha = %.2e, ', ...
                        'cut_lb = %.4e, nu = %.4e.\n'], ...
                        alpha_record(i),cut_lb_record(i),nu);
                end
                continue
            end
            rethrow(cut_error)
        end
        konno_record{i}.sdp_skipped=true;
        konno_record{i}.sdp_error_identifier=cut_error.identifier;
        konno_record{i}.sdp_error_message=cut_error.message;
        if parameters.verbose
            fprintf('  SDP cut failed; keeping the validated Konno cut and refreshing the relaxation next iteration.\n');
        end
        continue
    end

end

best_sol=best_sol*solution_scale;
if i>0
    bestlb=min([bestlb;cut_lb_record(1:i);konno_lb_record(1:i)]);
end
nb_iters=i;
bestub=bestub/parameters.scaling;
bestlb=bestlb/parameters.scaling;
lb_record=lb_record/parameters.scaling;
ub_record=ub_record/parameters.scaling;
bestub_record=bestub_record/parameters.scaling;
bestlb_record=bestlb_record/parameters.scaling;
cut_lb_record=cut_lb_record/parameters.scaling;
konno_lb_record=konno_lb_record(1:i)/parameters.scaling;
konno_record=konno_record(1:i); % Detailed records retain internal solver units.
time_record_konno=time_record_konno(1:i);
konno_recovery_record=konno_recovery_record(1:i);
time_record_konno_recovery=time_record_konno_recovery(1:i);
dc_recovery_record=dc_recovery_record(1:i);
time_record_dc_recovery=time_record_dc_recovery(1:i);
empty_retained_record=empty_retained_record(1:i);

lb_record=lb_record(1:i);
ub_record=ub_record(1:i);
bestub_record=bestub_record(1:i);
bestlb_record=bestlb_record(1:i);
cut_lb_record=cut_lb_record(1:i);
alpha_record=alpha_record(1:i);
time_record_lb=time_record_lb(1:i);
time_record_cut_lb=time_record_cut_lb(1:i);
time_record_kkt=time_record_kkt(1:i);
time_record_generate_cut=time_record_generate_cut(1:i);

% alpha_record is positive exactly when a DNN cut was appended. A
% generalized Konno record represents an appended cut only when it was
% accepted and its retained region was nonempty.
dnn_cut_count=nnz(alpha_record>0);
konno_cut_count=0;
for record_index=1:numel(konno_record)
    record=konno_record{record_index};
    if ~isempty(record) && isfield(record,'accepted') && record.accepted && ...
            isfield(record,'region_empty') && ~record.region_empty
        konno_cut_count=konno_cut_count+1;
    end
end
added_cut_count=size(A_bar,1)-initial_inequality_count;
recorded_cut_count=dnn_cut_count+konno_cut_count;
if recorded_cut_count~=added_cut_count
    warning('qpsolver:CutCountMismatch', ...
        ['Cut records report %d added cuts, but the retained system ', ...
         'contains %d added inequalities.'],recorded_cut_count,added_cut_count);
end
cut_counts=struct('total',added_cut_count,'dnn',dnn_cut_count, ...
    'konno',konno_cut_count,'records_match',recorded_cut_count==added_cut_count);
diagnostics=struct( ...
    'lb_record',lb_record, ...
    'ub_record',ub_record, ...
    'bestlb_record',bestlb_record, ...
    'bestub_record',bestub_record, ...
    'dnn_cut_lb_record',cut_lb_record, ...
    'alpha_record',alpha_record, ...
    'konno_record',{konno_record}, ...
    'konno_cut_lb_record',konno_lb_record, ...
    'konno_recovery_record',{konno_recovery_record}, ...
    'dc_recovery_record',{dc_recovery_record}, ...
    'empty_retained_record',{empty_retained_record}, ...
    'time_record_lb',time_record_lb, ...
    'time_record_ub',time_record_kkt, ...
    'time_record_generate_cut',time_record_generate_cut, ...
    'time_record_cut_verification',time_record_cut_lb, ...
    'time_record_konno',time_record_konno, ...
    'time_record_konno_recovery',time_record_konno_recovery, ...
    'time_record_dc_recovery',time_record_dc_recovery, ...
    'cut_counts',cut_counts, ...
    'konnofirst',konnofirst);

fprintf("************************************ End of Computation  ****************************\n");




end
