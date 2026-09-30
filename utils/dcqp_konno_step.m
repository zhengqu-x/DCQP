function [A_next,b_next,record] = dcqp_konno_step(Q,d,A,b,Aeq,beq,barx,nuR,r0,params)
%DCQP_KONNO_STEP Attempt a Konno cut at the current point without KKT tests.
% Input/output inequalities use DCQP's internal scaled coordinates.
A_next=A; b_next=b;
record=struct('attempted',true,'accepted',false,'region_empty',false, ...
    'status','not_started','message','','kkt_checked',false,'barx',barx, ...
    'c',[],'norm_c',NaN,'factor',NaN,'lower_bound',Inf,'target',nuR, ...
    'construction_status','','history',struct([]));
try
    % Paired inequalities let construction/validation retain any equalities
    % without a variable transformation. The main SDP still uses Aeq,beq.
    Ac=[A;Aeq;-Aeq]; bc=[b;beq;-beq];
    opts=struct('sdp_tolerance',params.mosek_tolerance, ...
        'verbose',logical(params.verbose),'mosek_quiet',params.mosek_quiet, ...
        'konno_options',struct('check_kkt',false,'verbose',logical(params.verbose), ...
        'lp_tolerance',params.gurobi_lp_tolerance));
    [c,info]=ValidatedExtendedKonnoCut(Q,d,Ac,bc,barx,nuR,r0,opts);
    record.status=info.status; record.message=info.message;
    record.construction_status=info.konno.status; record.history=info.history;
    if ~info.success || isempty(c), return; end
    record.c=c; record.norm_c=norm(c); record.factor=info.best_factor;
    if record.norm_c==0
        record.accepted=true; record.region_empty=true;
        record.lower_bound=info.best_lower_bound;
        record.status='whole_region_validated';
        return
    end
    candidate_A=[A;-c'/record.norm_c];
    candidate_b=[b;(-1-c'*barx)/record.norm_c];
    [~,~,flag]=gurobilp(zeros(numel(barx),1),candidate_A,candidate_b, ...
        Aeq,beq,[],[],params.gurobi_lp_method,params.gurobi_lp_tolerance);
    if flag==-1
        record.accepted=true; record.region_empty=true;
        record.lower_bound=info.best_lower_bound;
        record.status='retained_region_empty';
        return
    elseif flag~=1
        record.status='retained_region_check_failed';
        record.message='Could not check the retained region; skip the Konno cut.';
        return
    end
    [candidate_A(end,:),candidate_b(end),rescale_status]=rescale_constraint_by_slack( ...
        candidate_A(end,:),candidate_b(end),candidate_A,candidate_b,Aeq,beq, ...
        params.gurobi_lp_method,params.gurobi_lp_tolerance);
    if rescale_status==-1
        record.accepted=true; record.region_empty=true;
        record.lower_bound=info.best_lower_bound;
        record.status='retained_region_empty';
        return
    end
    A_next=candidate_A; b_next=candidate_b;
    record.accepted=true; record.lower_bound=info.best_lower_bound;
catch exception
    record.status='error'; record.message=exception.message;
end
end
