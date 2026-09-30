function [A,b,barx,v,bestub,best_sol,nu,nuR,beta,out] = ...
    compute_konno_cut(Q,d,A,b,Aeq,beq,M,N,barx,v,bestub,best_sol, ...
    nu,nuR,beta,barz,tstar,p)
%COMPUTE_KONNO_CUT Reproduce the stored proactive Konno attempt and retry.

out=struct('added',false,'refined',false,'record',[], ...
    'recovery',[],'recovery_time',0);
if p.verbose
    threshold=max(1e-5,abs(bestub)*p.gap_tolerance*0.1);
    fprintf(['  Konno first at current point (no KKT test), ', ...
        '|bestub-v|=%.4e <= %.2e\n'],abs(bestub-v),threshold);
end

[A,b,out.record]=dcqp_konno_step( ...
    Q,d,A,b,Aeq,beq,barx,nuR,abs(tstar),p);
out.record.relaxation_point_before=barz;
if ~out.record.accepted
    initial_attempt=out.record;
    if p.verbose
        fprintf('  Konno initial attempt failed: %s; %s\n', ...
            initial_attempt.status,initial_attempt.message);
    end
    recovery_timer=tic;
    [recovered_point,out.recovery]=refine_failed_cut_point( ...
        Q,d,A,b,Aeq,beq,M,N,barx,beta,p,'Konno',true);
    out.recovery_time=toc(recovery_timer);
    out.recovery.initial_attempt=initial_attempt;
    out.recovery.sdp_objective_point_before=barz;
    out.recovery.sdp_objective_point_after=barz;
    out.recovery.retry_attempted=false;
    out.recovery.retry_success=false;
    if out.recovery.refined
        out.refined=true;
        barx=recovered_point;
        v=barx'*Q*barx+2*d'*barx;
        if v<bestub
            best_sol=barx;
        end
        bestub=min(bestub,v);
        [nu,nuR,beta]=dcqp_cut_targets( ...
            bestub,v,p.eta,p.gap_tolerance);
        out.recovery.beta_after=beta;
        out.recovery.nu_after=nu;
        out.recovery.nuR_after=nuR;
        out.recovery.sdp_objective_point_after=barz;
        out.recovery.retry_attempted=true;
        if p.verbose
            fprintf(['  Konno recovery retry: bestub=%.4e, v=%.4e, ', ...
                'nu=%.4e, nu_R=%.4e, beta=%.2e; barz unchanged.\n'], ...
                bestub,v,nu,nuR,beta);
        end
        [A,b,out.record]=dcqp_konno_step( ...
            Q,d,A,b,Aeq,beq,barx,nuR,abs(tstar),p);
        out.record.relaxation_point_before=barz;
        out.recovery.retry_success=out.record.accepted;
    end
end

out.added=out.record.accepted && ~out.record.region_empty;
end
