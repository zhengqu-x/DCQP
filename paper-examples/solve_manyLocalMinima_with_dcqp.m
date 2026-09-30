function [summary,result_folder] = solve_manyLocalMinima_with_dcqp(selection,options)
%SOLVE_MANYLOCALMINIMA_WITH_DCQP Solve the many-local-minima StQP examples.
if nargin < 1, selection = "all"; end
if nargin < 2, options = struct(); end
[summary,result_folder] = run_paper_examples('structured','many_local_minima','dcqp',selection,options);
end
