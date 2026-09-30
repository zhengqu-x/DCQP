function [summary,result_folder] = solve_existing_testsets_with_gurobi(group_name,selection,options)
%SOLVE_EXISTING_TESTSETS_WITH_GUROBI Solve an existing-test group.
if nargin<2 || isempty(selection), selection="all"; end
if nargin<3, options=struct(); end
if isstruct(selection), options=selection; selection="all"; end
if nargin==2 && isnumeric(selection) && isscalar(selection) && selection>16
    options=struct('max_time',selection); selection="all"; % Legacy Tlimit.
elseif isnumeric(options) && isscalar(options)
    options=struct('max_time',options);
end
[summary,result_folder]=run_paper_examples( ...
    'existing-tests',group_name,'gurobi',selection,options);
end
