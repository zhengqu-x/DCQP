function [summary,result_folder] = solve_boxQP_with_dcqp(selection,options)
%SOLVE_BOXQP_WITH_DCQP Solve the structured BoxQP examples.
if nargin < 1, selection = "all"; end
if nargin < 2, options = struct(); end
[summary,result_folder] = run_paper_examples('structured','boxqp','dcqp',selection,options);
end
