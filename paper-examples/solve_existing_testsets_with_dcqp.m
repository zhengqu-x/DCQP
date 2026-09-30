function [summary,result_folder] = solve_existing_testsets_with_dcqp(group_name,selection,options)
%SOLVE_EXISTING_TESTSETS_WITH_DCQP Solve an existing-test group.
if nargin<2 || isempty(selection), selection="all"; end
if nargin<3, options=struct(); end
if isstruct(selection), options=selection; selection="all"; end
[summary,result_folder]=run_paper_examples( ...
    'existing-tests',group_name,'dcqp',selection,options);
end
