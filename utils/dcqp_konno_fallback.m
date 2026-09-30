function [A,b,state,record,recovery] = dcqp_konno_fallback(Q,d,A,b,Aeq,beq,M,N,state,barz,r0,p,check_residual,allow_retry)
%DCQP_KONNO_FALLBACK Try a validated Konno cut after an SDP cut failure.
% The caller retains the SDP relaxation point barz and handles global bounds.
% check_residual is false when SDP failure recovery already checked/refined
% the current anchor. allow_retry prevents a second Konno recovery batch in
% an outer iteration that already performed its unconditional 100 steps.
recovery=[]; residual_check=[];
if p.verbose
    fprintf('SDP failure: activating Konno without the konnofirst/closeness gate.\n');
end
if check_residual
    [point,residual_check]=refine_failed_cut_point(Q,d,A,b,Aeq,beq, ...
        M,N,state.barx,state.beta,p);
    if residual_check.refined
        state=refresh_state(state,point,Q,d,p);
        recovery=residual_check;
        recovery.sdp_objective_point_before=barz;
        recovery.sdp_objective_point_after=barz;
    end
end
[A,b,record]=dcqp_konno_step(Q,d,A,b,Aeq,beq,state.barx,state.nuR,r0,p);
if ~record.accepted && allow_retry
    first_attempt=record;
    [point,konno_recovery]=refine_failed_cut_point(Q,d,A,b,Aeq,beq, ...
        M,N,state.barx,state.beta,p,'Konno',true);
    konno_recovery.initial_attempt=first_attempt;
    konno_recovery.sdp_objective_point_before=barz;
    konno_recovery.sdp_objective_point_after=barz;
    konno_recovery.retry_attempted=false;
    konno_recovery.retry_success=false;
    konno_recovery.initial_residual_check=residual_check;
    if konno_recovery.refined
        state=refresh_state(state,point,Q,d,p);
        konno_recovery.retry_attempted=true;
        [A,b,record]=dcqp_konno_step(Q,d,A,b,Aeq,beq,state.barx,state.nuR,r0,p);
        konno_recovery.retry_success=record.accepted;
    end
    recovery=konno_recovery;
end
record.sdp_failure_fallback=true;
record.kkt_residual_check=residual_check;
record.relaxation_point_before=barz;
record.sdp_objective_point=barz;
record.sdp_region_rows=size(A,1);
if ~record.accepted && p.verbose
    fprintf('Fallback Konno skipped: %s; %s\n',record.status,record.message);
end
end

function state=refresh_state(state,point,Q,d,p)
state.barx=point;
state.v=point'*Q*point+2*d'*point;
if state.v<state.bestub, state.best_sol=point; end
state.bestub=min(state.bestub,state.v);
[state.nu,state.nuR,state.beta]=dcqp_cut_targets( ...
    state.bestub,state.v,p.eta,p.gap_tolerance);
end
