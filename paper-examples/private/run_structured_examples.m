function [summary,result_folder] = run_structured_examples(family,solver,selection,options)
%RUN_STRUCTURED_EXAMPLES Compatibility wrapper for the unified runner.
if nargin<3 || isempty(selection), selection="all"; end
if nargin<4, options=struct(); end
[summary,result_folder]=run_paper_examples( ...
    'structured',family,solver,selection,options);
end
