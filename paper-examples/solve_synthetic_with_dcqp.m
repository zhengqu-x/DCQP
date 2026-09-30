function [summary,result_folder] = solve_synthetic_with_dcqp(group_name,selection,options)
%SOLVE_SYNTHETIC_WITH_DCQP Solve a synthetic group using unified results.
if nargin<2 || isempty(selection), selection="all"; end
if nargin<3, options=struct(); end
if isstruct(selection), options=selection; selection="all"; end
[summary,result_folder]=run_paper_examples( ...
    'synthetic',group_name,'dcqp',selection,options);
end
