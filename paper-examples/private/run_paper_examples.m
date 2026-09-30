function [summary,result_folder] = run_paper_examples(category,family,solver,selection,options)
%RUN_PAPER_EXAMPLES Shared runner and result schema for all paper examples.

if nargin < 4 || isempty(selection), selection = "all"; end
if nargin < 5, options = struct(); end
validateattributes(category,{'char','string'},{'scalartext'});
validateattributes(family,{'char','string'},{'scalartext'});
validateattributes(solver,{'char','string'},{'scalartext'});
category=char(category); family=char(family); solver=char(solver);

categories={'structured','synthetic','existing-tests'};
solvers={'dcqp','dcqp_eliminate','gurobi'};
assert(ismember(category,categories),'paper_examples:category', ...
    'Unknown experiment category: %s',category);
assert(ismember(solver,solvers),'paper_examples:solver', ...
    'Unknown solver: %s',solver);
assert(~strcmp(solver,'dcqp_eliminate') || strcmp(category,'synthetic'), ...
    'paper_examples:solver','dcqp_eliminate is only available for synthetic instances.');
assert(isstruct(options) && isscalar(options),'paper_examples:options', ...
    'options must be a scalar structure.');

[data_folder,pattern,defaults]=configuration(category,family);
unknown=setdiff(fieldnames(options),fieldnames(defaults));
if ~isempty(unknown)
    error('paper_examples:options','Unknown option: %s',unknown{1});
end
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(options,names{k}), options.(names{k})=defaults.(names{k}); end
end

files=dir(fullfile(data_folder,pattern));
files=natural_sort(files);
files=select_files(files,selection);
assert(~isempty(files),'paper_examples:selection','No instances were selected.');

paper_folder=fileparts(fileparts(mfilename('fullpath')));
root=fileparts(paper_folder);
addpath(root,fullfile(root,'utils'),fullfile(root,'legacy'));
result_folder=fullfile(paper_folder,'results',category,family, ...
    sprintf('%s_%s',solver,char(datetime('now','Format','yyyyMMdd_HHmmss'))));
mkdir(result_folder);
run_started=datetime('now','Format','yyyy-MM-dd HH:mm:ss');

if options.write_log
    diary(fullfile(result_folder,'run.log'));
    diary on;
    diary_cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
end

blank=struct('id','','category',category,'family',family,'solver',solver, ...
    'original_n',NaN,'seed',NaN,'status','pending','solved',false, ...
    'time',NaN,'initial_relative_gap',NaN,'relative_gap',Inf, ...
    'iterations',NaN,'added_cuts',NaN,'upper_bound',Inf, ...
    'lower_bound',-Inf,'feasibility_violation',Inf,'error_message','');
rows=repmat(blank,numel(files),1);

for k=1:numel(files)
    id=erase(files(k).name,'.mat');
    [instance,meta]=load_instance(fullfile(files(k).folder,files(k).name), ...
        category,family,id);
    row=blank; row.id=id;
    row.original_n=get_field(meta,'original_n',size(instance.Q,1));
    row.seed=get_field(meta,'seed',NaN);
    fprintf('[%d/%d] %s with %s\n',k,numel(files),id,upper(solver));
    timer=tic; x=[]; raw=struct();
    try
        if isfinite(row.seed), rng(row.seed,'twister'); end
        if startsWith(solver,'dcqp')
            params=dcqp_parameters(options,id,meta);
            if strcmp(solver,'dcqp_eliminate') && ~isempty(instance.Aeq)
                n=size(instance.Q,1);
                [Q2,d2,A2,b2,M,x0,equality_constant]=eliminate_equalities( ...
                    instance.Q,instance.c/2,instance.A,instance.b, ...
                    instance.Aeq,instance.beq,n);
                [y,~,raw]=dcqp_solve(Q2,d2,A2,b2,[],[],params);
                x=M*y+x0;
                raw=restore_objective_constant(raw,equality_constant);
            else
                [x,~,raw]=dcqp_solve(instance.Q,instance.c/2, ...
                    instance.A,instance.b,instance.Aeq,instance.beq,params);
            end
            if strcmp(raw.status,'not solved') && isfinite(options.retry_eta)
                first_attempt=raw;
                params.eta=options.retry_eta;
                if strcmp(solver,'dcqp_eliminate') && ~isempty(instance.Aeq)
                    [y,~,raw]=dcqp_solve(Q2,d2,A2,b2,[],[],params);
                    x=M*y+x0;
                    raw=restore_objective_constant(raw,equality_constant);
                else
                    [x,~,raw]=dcqp_solve(instance.Q,instance.c/2, ...
                        instance.A,instance.b,instance.Aeq,instance.beq,params);
                end
                raw.first_attempt=first_attempt;
            end
            row.status=raw.status;
            row.initial_relative_gap=raw.initial_relative_gap;
            row.relative_gap=raw.gap;
            row.iterations=raw.iterations;
            row.added_cuts=raw.added_cuts;
            row.upper_bound=raw.upper_bound;
            row.lower_bound=raw.lower_bound;
            if isfield(raw,'error_message'), row.error_message=raw.error_message; end
        else
            [x,raw]=solve_with_gurobi(instance,meta,options);
            row.status=raw.status;
            if isfield(raw,'objval'), row.upper_bound=raw.objval; end
            if isfield(raw,'objbound'), row.lower_bound=raw.objbound; end
            if isfield(raw,'mipgap')
                row.relative_gap=raw.mipgap;
            elseif isfinite(row.upper_bound) && isfinite(row.lower_bound)
                row.relative_gap=relative_gap(row.upper_bound,row.lower_bound, ...
                    options.gap_tolerance);
            end
        end
    catch ME
        row.status='error'; row.error_message=ME.message;
        raw=struct('status','error','error_message',ME.message, ...
            'error_stack',ME.stack);
    end
    row.time=toc(timer);
    if numel(x)==size(instance.Q,1) && all(isfinite(x))
        row.feasibility_violation=max([0;instance.A*x-instance.b; ...
            abs(instance.Aeq*x-instance.beq);instance.LB-x;x-instance.UB]);
        row.upper_bound=x'*instance.Q*x+instance.c'*x;
    end
    row.solved=row.relative_gap<=options.gap_tolerance && ...
        row.relative_gap>=-1e-8 && ...
        row.feasibility_violation<=options.feasibility_tolerance && ...
        ~strcmp(row.status,'error');
    rows(k)=row;
    save(fullfile(result_folder,[id '.mat']),'row','raw','x','options','meta');
    try
        write_live_summary(result_folder,rows(1:k),numel(files),run_started);
    catch summary_error
        warning('paper_examples:liveSummary', ...
            'Could not update summary.md: %s',summary_error.message);
    end
end

summary=struct2table(rows,'AsArray',true);
writetable(summary,fullfile(result_folder,'summary.csv'));
save(fullfile(result_folder,'summary.mat'),'summary','rows','options', ...
    'category','family','solver');
fprintf('Finished %d cases; %d reached the target. Results: %s\n', ...
    numel(rows),sum([rows.solved]),result_folder);
end

function [data_folder,pattern,defaults]=configuration(category,family)
paper_folder=fileparts(fileparts(mfilename('fullpath')));
root=fileparts(paper_folder);
defaults=struct('max_time',3600,'gap_tolerance',1e-4, ...
    'max_iterations',300,'nb_rounds',1,'do_scaling',true, ...
    'konnofirst',false,'known_solution_start',false,'verbose',true, ...
    'mosek_quiet',true,'mosek_tolerance',1e-8, ...
    'accept_cut_below_threshold',false,'retry_eta',NaN, ...
    'gurobi_threads',0,'feasibility_tolerance',1e-6,'write_log',true);
switch category
    case 'structured'
        valid={'inexact_stqp','boxqp','many_local_minima'};
        assert(ismember(family,valid),'paper_examples:family', ...
            'Unknown structured family: %s',family);
        data_folder=fullfile(root,'data','structured',family);
        pattern='*.mat';
        defaults.max_time=600;
        defaults.konnofirst=true;
        defaults.known_solution_start=true;
    case 'synthetic'
        valid={'qp_n_0_1','qp_n_0_3','qp_n_0_9', ...
            'qp_u_0_1','qp_u_0_3','qp_u_0_9','qp_u_25_1'};
        assert(ismember(family,valid),'paper_examples:family', ...
            'Unknown synthetic family: %s',family);
        data_folder=fullfile(root,'data','synthetic');
        pattern=[family '_*.mat'];
        defaults.nb_rounds=10;
        defaults.mosek_tolerance=1e-9;
        defaults.retry_eta=0.8;
    case 'existing-tests'
        valid={'qp20_10','qp30_15','qp40_20','qp50_25'};
        assert(ismember(family,valid),'paper_examples:family', ...
            'Unknown existing-test family: %s',family);
        data_folder=fullfile(root,'data','existing_testsets');
        pattern=[family '_*.mat'];
end
assert(isfolder(data_folder),'paper_examples:data', ...
    'Data folder not found: %s',data_folder);
end

function params=dcqp_parameters(options,id,meta)
params=dcqp_default_params();
params.max_time=options.max_time;
params.gap_tolerance=options.gap_tolerance;
params.max_iterations=options.max_iterations;
params.nb_rounds=options.nb_rounds;
params.do_scaling=options.do_scaling;
params.konnofirst=options.konnofirst;
params.mosek_tolerance=options.mosek_tolerance;
params.accept_cut_below_threshold=options.accept_cut_below_threshold;
params.verbose=options.verbose;
params.mosek_quiet=options.mosek_quiet;
params.display_summary=options.verbose;
params.robust_mode=true;
params.save_failed_instances=false;
params.filename=id;
if options.known_solution_start && isfield(meta,'known_solution') && ...
        ~isempty(meta.known_solution)
    params.known_solution=meta.known_solution;
end
end

function [x,raw]=solve_with_gurobi(instance,meta,options)
model=struct('Q',sparse(instance.Q),'obj',full(instance.c(:)), ...
    'A',sparse([instance.A;instance.Aeq]), ...
    'rhs',full([instance.b(:);instance.beq(:)]), ...
    'sense',[repmat('<',size(instance.A,1),1); ...
             repmat('=',size(instance.Aeq,1),1)], ...
    'lb',-Inf(size(instance.Q,1),1),'ub',Inf(size(instance.Q,1),1), ...
    'modelsense','min');
if options.known_solution_start && isfield(meta,'known_solution') && ...
        ~isempty(meta.known_solution)
    model.start=full(meta.known_solution(:));
end
% Match the numerical settings used by the original paper-example scripts.
gp=struct('OptimalityTol',options.gap_tolerance, ...
    'BarConvTol',options.gap_tolerance,'LPWarmStart',2, ...
    'TimeLimit',options.max_time,'OutputFlag',double(options.verbose));
if options.gurobi_threads>0, gp.Threads=options.gurobi_threads; end
raw=gurobi(model,gp);
if strcmp(raw.status,'INF_OR_UNBD')
    gp.DualReductions=0;
    raw=gurobi(model,gp);
end
if isfield(raw,'x'), x=raw.x; else, x=[]; end
end

function [instance,meta]=load_instance(file,category,family,id)
loaded=load(file);
if strcmp(category,'structured')
    instance=loaded.instance;
    if isfield(instance,'metadata'), meta=instance.metadata;
    else, meta=struct(); end
else
    n=size(loaded.H,1);
    Aeq=loaded.Aeq; beq=loaded.beq;
    if isempty(Aeq) || norm(Aeq,'fro')==0, Aeq=[]; beq=[]; end
    LB=full(loaded.LB(:)); UB=full(loaded.UB(:));
    instance=struct('Q',loaded.H/2,'c',full(loaded.f(:)), ...
        'A',full([loaded.A;eye(n);-eye(n)]), ...
        'b',full([loaded.b(:);UB;-LB]),'Aeq',full(Aeq), ...
        'beq',full(beq(:)),'LB',LB,'UB',UB);
    meta=struct('family',family,'source_id',id,'original_n',n);
end
required={'Q','c','A','b','Aeq','beq','LB','UB'};
for k=1:numel(required)
    assert(isfield(instance,required{k}),'paper_examples:data', ...
        'Instance %s is missing field %s.',id,required{k});
end
instance.Q=full(instance.Q); instance.c=full(instance.c(:));
instance.A=full(instance.A); instance.b=full(instance.b(:));
if isempty(instance.Aeq)
    instance.Aeq=zeros(0,size(instance.Q,1));
    instance.beq=zeros(0,1);
else
    instance.Aeq=full(instance.Aeq);
    instance.beq=full(instance.beq(:));
end
instance.LB=full(instance.LB(:)); instance.UB=full(instance.UB(:));
end

function raw=restore_objective_constant(raw,constant)
fields={'upper_bound','lower_bound','initial_upper_bound','initial_lower_bound'};
for k=1:numel(fields)
    if isfield(raw,fields{k}) && isfinite(raw.(fields{k}))
        raw.(fields{k})=raw.(fields{k})+constant;
    end
end
raw.equality_elimination_constant=constant;
end

function value=get_field(s,name,default)
if isfield(s,name) && ~isempty(s.(name)), value=s.(name);
else, value=default; end
end

function gap=relative_gap(upper,lower,tolerance)
if abs(upper)>tolerance, gap=abs(upper-lower)/abs(upper);
else, gap=abs(upper-lower); end
end

function files=natural_sort(files)
if isempty(files), return; end
keys=cell(size(files));
for k=1:numel(files)
    parts=regexp(files(k).name,'\d+|\D+','match');
    key='';
    for j=1:numel(parts)
        if ~isempty(regexp(parts{j},'^\d+$','once'))
            key=[key sprintf('%012d',str2double(parts{j}))]; %#ok<AGROW>
        else
            key=[key lower(parts{j})]; %#ok<AGROW>
        end
    end
    keys{k}=key;
end
[~,order]=sort(keys); files=files(order);
end

function selected=select_files(files,selection)
if isstring(selection) && isscalar(selection) && strcmpi(selection,"all")
    selected=files; return
end
if ischar(selection) && strcmpi(selection,'all')
    selected=files; return
end
if isnumeric(selection)
    validateattributes(selection,{'numeric'},{'vector','integer','positive'});
    assert(all(selection<=numel(files)),'paper_examples:selection', ...
        'Selection index exceeds the number of available instances.');
    selected=files(selection); return
end
requested=cellstr(string(selection));
ids=erase({files.name},'.mat');
[found,positions]=ismember(requested,ids);
assert(all(found),'paper_examples:selection','Unknown instance: %s', ...
    requested{find(~found,1)});
selected=files(positions);
end

function write_live_summary(result_folder,rows,expected_count,started)
completed=numel(rows);
solved=sum([rows.solved]);
errors=sum(strcmp({rows.status},'error'));
target=fullfile(result_folder,'summary.md');
temporary=[target '.tmp'];
fid=fopen(temporary,'w');
assert(fid>=0,'paper_examples:liveSummary','Cannot write %s.',temporary);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'# Live experiment summary\n\n');
fprintf(fid,'- Result directory: `%s`\n',result_folder);
fprintf(fid,'- Started: %s\n',char(started));
refreshed=datetime('now','Format','yyyy-MM-dd HH:mm:ss');
fprintf(fid,'- Last refresh: %s\n',char(refreshed));
fprintf(fid,'- Completed: **%d / %d**\n',completed,expected_count);
fprintf(fid,'- Reached target: **%d**\n',solved);
fprintf(fid,'- Errors: **%d**\n\n',errors);
fprintf(fid,['| # | Instance | Status | Solved | Time (s) | Initial gap | ', ...
    'Relative gap | Iter. | Added cuts | Upper bound | Lower bound | Feasibility |\n']);
fprintf(fid,'|---:|---|---|:---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for k=1:completed
    row=rows(k);
    fprintf(fid,['| %d | %s | %s | %s | %.3f | %.3e | %.3e | %.0f | ', ...
        '%.0f | %.8e | %.8e | %.3e |\n'], ...
        k,md_escape(row.id),md_escape(row.status),yes_no(row.solved), ...
        row.time,row.initial_relative_gap,row.relative_gap,row.iterations, ...
        row.added_cuts,row.upper_bound,row.lower_bound, ...
        row.feasibility_violation);
end
clear cleanup
movefile(temporary,target,'f');
end

function value=md_escape(value)
value=char(string(value));
value=strrep(value,'|','\|');
value=strrep(value,newline,' ');
end

function value=yes_no(flag)
if flag, value='yes'; else, value='no'; end
end
