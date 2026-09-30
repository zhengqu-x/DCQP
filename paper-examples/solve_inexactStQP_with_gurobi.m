function [summary,result_folder] = solve_inexactStQP_with_gurobi(selection,options)
%SOLVE_INEXACTSTQP_WITH_GUROBI Solve the structured inexact StQP examples.
if nargin < 1, selection = "all"; end
if nargin < 2, options = struct(); end
[summary,result_folder] = run_paper_examples('structured','inexact_stqp','gurobi',selection,options);
end
