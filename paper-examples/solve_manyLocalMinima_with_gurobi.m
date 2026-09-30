function [summary,result_folder] = solve_manyLocalMinima_with_gurobi(selection,options)
%SOLVE_MANYLOCALMINIMA_WITH_GUROBI Solve the many-local-minima StQP examples.
if nargin < 1, selection = "all"; end
if nargin < 2, options = struct(); end
[summary,result_folder] = run_paper_examples('structured','many_local_minima','gurobi',selection,options);
end
