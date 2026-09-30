function [summary,result_folder] = solve_boxQP_with_gurobi(selection,options)
%SOLVE_BOXQP_WITH_GUROBI Solve the structured BoxQP examples.
if nargin < 1, selection = "all"; end
if nargin < 2, options = struct(); end
[summary,result_folder] = run_paper_examples('structured','boxqp','gurobi',selection,options);
end
