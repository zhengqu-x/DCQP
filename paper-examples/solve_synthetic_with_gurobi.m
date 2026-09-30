function [summary,result_folder] = solve_synthetic_with_gurobi(group_name,selection,options)
%SOLVE_SYNTHETIC_WITH_GUROBI Solve a synthetic group using unified results.
if nargin<2 || isempty(selection), selection="all"; end
if nargin<3, options=struct(); end
if isstruct(selection), options=selection; selection="all"; end
if isnumeric(options) && isscalar(options)
    options=struct('max_time',options); % Backward-compatible Tlimit argument.
end
[summary,result_folder]=run_paper_examples( ...
    'synthetic',group_name,'gurobi',selection,options);
end
