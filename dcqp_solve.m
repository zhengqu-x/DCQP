function [x_opt, fval, info] = dcqp_solve(Q, d, A, b, Aeq, beq, params)
% DCQP_SOLVE  Solve nonconvex quadratic programming problems using doubly nonnegative cutting planes
%
% SYNTAX:
%   [x_opt, fval, info] = dcqp_solve(Q, d, A, b)
%   [x_opt, fval, info] = dcqp_solve(Q, d, A, b, Aeq, beq)
%   [x_opt, fval, info] = dcqp_solve(Q, d, A, b, Aeq, beq, params)
%
% DESCRIPTION:
%   Solves the quadratic programming problem:
%       minimize    x'*Q*x + 2*d'*x
%       subject to  A*x <= b
%                   Aeq*x = beq
%
%
%   REQUIREMENTS:
%   - The matrix Q has at least one negative eigenvalue
%   - Inequality constraints A*x <= b are MANDATORY (cannot be empty)
%   - The feasible region is BOUNDED
%
% INPUT:
%   Q       - (n x n) symmetric matrix (objective Hessian)
%   d       - (n x 1) vector (linear objective coefficients)
%   A       - (m x n) matrix (inequality constraint coefficients, REQUIRED)
%   b       - (m x 1) vector (inequality constraint bounds, REQUIRED)
%   Aeq     - (meq x n) matrix (equality constraint coefficients, optional)
%   beq     - (meq x 1) vector (equality constraint bounds, optional)
%   params  - Structure with algorithm parameters (optional)
%             .known_solution: known optimal point in input coordinates;
%             when supplied, initialize directly without random starts.
%
% OUTPUT:
%   x_opt   - (n x 1) optimal solution vector
%   fval    - Optimal objective function value
%   info    - Structure with solution information:
%             .status      - Solution status string
%             .gap         - Relative optimality gap
%             .iterations  - Number of iterations
%             .time        - Total computation time
%             .upper_bound - Best upper bound found
%             .lower_bound - Best lower bound found
%             .added_cuts  - Number of inequalities added by DCQP
%             .dnn_cuts    - Number of added DNN cuts
%             .konno_cuts  - Number of added generalized Konno cuts
%             .initial_relative_gap - Relative gap after the first bound
%             .diagnostics - Detailed per-iteration records
%
% EXAMPLES:
%   % Simple nonconvex QP
%   Q = [1 -1; -1 -1]; d = [1; 1]; A = [1 1]; b = 1;
%   [x, fval, info] = dcqp_solve(Q, d, A, b);
%
%   % With equality constraints
%   Aeq = [1 0]; beq = 0.5;
%   [x, fval, info] = dcqp_solve(Q, d, A, b, Aeq, beq);
%
%   % With custom parameters
%   params = dcqp_default_params();
%   params.gap_tolerance = 1e-6;
%   [x, fval, info] = dcqp_solve(Q, d, A, b, [], [], params);
%
% SEE ALSO: dcqp_default_params, qpsolver, DC_decomposition

% Copyright (c) 2025
% All rights reserved.

% Input validation and preprocessing
if nargin < 4
    error('dcqp_solve:insufficient_inputs', ...
          'At least Q, d, A, and b must be provided');
end

% Set default values for optional inputs
if nargin < 5 || isempty(Aeq)
    Aeq = [];
end
if nargin < 6 || isempty(beq)
    beq = [];
end
if nargin < 7 || isempty(params)
    params = dcqp_default_params();
end





try
    % Validate inputs
    maintimer=tic;
    [n,Q,d,A,b,Aeq,beq,x_shift,obj_constant]=dcqp_check_input(Q, d, A, b, Aeq, beq, params);
    
    nb_rounds = max(1, min(n, params.nb_rounds));
    use_known_solution=isfield(params,'known_solution') && ~isempty(params.known_solution);
    if use_known_solution
        validateattributes(params.known_solution,{'double'},{'real','finite','vector','numel',n});
        sol=full(params.known_solution(:))-x_shift;
        violation=max([0;A*sol-b]);
        if ~isempty(Aeq), violation=max(violation,max(abs(Aeq*sol-beq))); end
        if violation>max(1e-8,10*params.gurobi_lp_tolerance)
            error('dcqp_solve:invalid_known_solution','The supplied known solution is infeasible (violation %.2e).',violation);
        end
        ub=sol'*Q*sol+2*d'*sol;
        initialization='known_solution';
    else
        [ub,sol]=compute_ub(Q,d,A,b,Aeq,beq,n,params,[],nb_rounds);
        initialization='random';
    end

    if params.do_scaling==true
        if abs(ub)>=5
            params.scaling=1/floor(abs(ub));
        elseif abs(ub)<=0.5 && abs(ub)>1e-3
            params.scaling=ceil(1/abs(ub));
        end

        Q=Q*params.scaling;
        d=d*params.scaling;
    end
    
   
    if use_known_solution
        sol2=sol; ub2=sol2'*Q*sol2+2*d'*sol2;
    else
        [ub2,sol2]=compute_ub(Q,d,A,b,Aeq,beq,n,params,sol,nb_rounds);
    end

    if params.verbose==true
        fprintf('initialization = %s\n',initialization);
        fprintf('initial objective value = %.4e\n',ub2/params.scaling+obj_constant);
        fprintf('initial bestub (solver units) = %.4e, scaling = %.2e, shift constant = %.2e\n', ...
            ub2, params.scaling, obj_constant);
    end
    params.solve_timer=maintimer;
    [best_ub,best_sol,best_lb,nb_iters,cut_counts,diagnostics]=qpsolver( ...
        Q,d,A,b,Aeq,beq,-Inf,ub2,sol2,params);
    total_time=toc(maintimer);
    
    x_opt = best_sol + x_shift;
    fval = best_ub + obj_constant;
    final_lb = best_lb + obj_constant;
    
    % Compute final gap in shifted solver coordinates to avoid sensitivity
    % to the additive objective constant introduced by variable shifting.
    if abs(best_ub) > params.gap_tolerance
        gap = abs(best_ub - best_lb) / abs(best_ub);
    else
        gap = abs(best_ub - best_lb);
    end
    if abs(fval) > params.gap_tolerance
        original_gap = abs(fval - final_lb) / abs(fval);
    else
        original_gap = abs(fval - final_lb);
    end

    if isempty(diagnostics.bestub_record) || isempty(diagnostics.bestlb_record)
        initial_upper_bound=best_ub;
        initial_lower_bound=best_lb;
    else
        initial_upper_bound=diagnostics.bestub_record(1);
        initial_lower_bound=diagnostics.bestlb_record(1);
    end
    if abs(initial_upper_bound)>params.gap_tolerance
        initial_gap=abs(initial_upper_bound-initial_lower_bound)/ ...
            abs(initial_upper_bound);
    else
        initial_gap=abs(initial_upper_bound-initial_lower_bound);
    end
    
    

    % Determine solution status
    if gap <= params.gap_tolerance
        status = ['successfully reduced relative gap below ' num2str(params.gap_tolerance)] ;
    elseif total_time >= params.max_time
        status = 'time_limit reached';
    elseif nb_iters>=params.max_iterations
        status = 'iteration_limit reached';
    else
        status='not solved';
    end
    
    % Package solution information
    info = struct();
    info.status = status;
    info.gap = gap;
    info.iterations = nb_iters; 
    info.time = total_time;
    info.upper_bound = fval;
    info.lower_bound = final_lb;
    info.initial_relative_gap = initial_gap;
    info.initial_upper_bound = initial_upper_bound + obj_constant;
    info.initial_lower_bound = initial_lower_bound + obj_constant;
    info.added_cuts = cut_counts.total;
    info.dnn_cuts = cut_counts.dnn;
    info.konno_cuts = cut_counts.konno;
    info.cut_records_match = cut_counts.records_match;
    info.original_gap = original_gap;
    info.scaling = params.scaling;
    info.variable_shift = x_shift;
    info.objective_constant = obj_constant;
    info.initialization = initialization;
    info.initial_solution = sol2+x_shift;
    info.diagnostics = diagnostics;

    
    
    % Display summary if verbose
    if params.display_summary
        fprintf('\n=== DCQP Solution Summary ===\n');
        fprintf('Instance name='); fprintf(params.filename);fprintf('\n');
        fprintf('Status: %s\n', status);
        fprintf('Best objective value: %.4e\n', fval);
        fprintf('Relative gap: %.2e\n', gap);
        fprintf('Original-coordinate relative gap: %.2e\n', original_gap);
        fprintf('Variable shift norm: %.2e\n', norm(x_shift));
        fprintf('Objective constant from shift: %.2e\n', obj_constant);
        fprintf('Computation time: %.2f seconds\n', total_time);
        fprintf('Number of iterations: %d\n', info.iterations);
        fprintf('Number of added cuts: %d (DNN: %d, generalized Konno: %d)\n', ...
            info.added_cuts,info.dnn_cuts,info.konno_cuts);
        fprintf('Problem dimension: %d variables, %d inequality constraints, %d equality constraints\n', n, size(A,1),size(Aeq,1));
        fprintf('==============================\n\n');
    end
    
catch ME
    % Handle solver errors gracefully
    total_time = toc(maintimer);

    x_opt = [];
    fval = inf;
    info = struct();
    info.status = 'error';
    info.gap = inf;
    info.iterations = 0;
    info.time = total_time;
    info.upper_bound = inf;
    info.lower_bound = -inf;
    info.initial_relative_gap = NaN;
    info.initial_upper_bound = Inf;
    info.initial_lower_bound = -Inf;
    info.added_cuts = NaN;
    info.dnn_cuts = NaN;
    info.konno_cuts = NaN;
    info.cut_records_match = false;
    info.diagnostics = struct();
    info.error_message = ME.message;
    info.error_stack = ME.stack;

    info

    % Save failed instance for debugging (use original inputs)
    if isfield(params, 'save_failed_instances') && params.save_failed_instances
        save_failed_instance(Q,d,A,b,Aeq,beq, params, ME, total_time);
    end

    if params.display_summary
        fprintf('\n=== DCQP Error ===\n');
        fprintf('Error: %s\n', ME.message);
        fprintf('==================\n\n');
    end

    % Re-throw error if not in robust mode
    if ~params.robust_mode
        rethrow(ME);
    end
end

end
