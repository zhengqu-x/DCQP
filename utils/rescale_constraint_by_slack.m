function [a_row,b_value] = rescale_constraint_by_slack(a_row,b_value,A,b,Aeq,beq,met_glp,tol_glp)

[~,min_row_value,exitflag]=gurobilp(a_row',A,b,Aeq,beq,[],[],met_glp,tol_glp);
if exitflag~=1
    error('Failed to compute constraint slack bound while rescaling constraints.');
end

slack_bound=b_value-min_row_value;
if slack_bound < -1e-8
    error('Computed a negative constraint slack bound while rescaling constraints.');
end
slack_bound=max(slack_bound,0);

if slack_bound>0.1
    scale=0.1/slack_bound;
    a_row=a_row*scale;
    b_value=b_value*scale;
end

end
