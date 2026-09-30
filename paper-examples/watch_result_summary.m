function watch_result_summary(result_folder,expected_count,poll_seconds)
%WATCH_RESULT_SUMMARY Continuously summarize any unified paper-example run.

if nargin<2 || isempty(expected_count), expected_count=number_of_cases(result_folder); end
if nargin<3 || isempty(poll_seconds), poll_seconds=10; end
validateattributes(expected_count,{'numeric'},{'scalar','integer','positive'});
validateattributes(poll_seconds,{'numeric'},{'scalar','positive'});
assert(isfolder(result_folder),'watch_result_summary:folder', ...
    'Result folder does not exist: %s',result_folder);

started=datetime('now','Format','yyyy-MM-dd HH:mm:ss');
while true
    completed=write_live_summary(result_folder,expected_count,started);
    if completed>=expected_count || isfile(fullfile(result_folder,'summary.mat'))
        write_live_summary(result_folder,expected_count,started);
        return
    end
    pause(poll_seconds);
end
end

function count=number_of_cases(result_folder)
summary_file=fullfile(result_folder,'summary.mat');
if isfile(summary_file)
    saved=load(summary_file,'rows'); count=numel(saved.rows); return
end
normalized=strrep(result_folder,'\','/');
if contains(normalized,'/synthetic/'), count=20;
elseif contains(normalized,'/existing-tests/'), count=16;
elseif contains(normalized,'/structured/many_local_minima/'), count=15;
else, count=120;
end
end

function completed=write_live_summary(result_folder,expected_count,started)
files=dir(fullfile(result_folder,'*.mat'));
files=files(~ismember({files.name},{'summary.mat'}));
files=natural_sort(files);
rows=struct([]);
for k=1:numel(files)
    try
        saved=load(fullfile(files(k).folder,files(k).name),'row');
        if isfield(saved,'row')
            if isempty(rows), rows=saved.row; else, rows(end+1)=saved.row; end %#ok<AGROW>
        end
    catch
        % Retry a MAT file on the next refresh if it is still being written.
    end
end
completed=numel(rows);
solved=0; errors=0;
if completed>0
    solved=sum([rows.solved]); errors=sum(strcmp({rows.status},'error'));
end

target=fullfile(result_folder,'summary.md'); temporary=[target '.tmp'];
fid=fopen(temporary,'w');
assert(fid>=0,'watch_result_summary:file','Cannot write %s.',temporary);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'# Live experiment summary\n\n');
fprintf(fid,'- Result directory: `%s`\n',result_folder);
fprintf(fid,'- Started watching: %s\n',char(started));
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
    initial_gap=get_field(row,'initial_relative_gap',NaN);
    final_gap=get_field(row,'relative_gap',get_field(row,'solver_gap',NaN));
    added_cuts=get_field(row,'added_cuts',NaN);
    fprintf(fid,['| %d | %s | %s | %s | %.3f | %.3e | %.3e | %.0f | ', ...
        '%.0f | %.8e | %.8e | %.3e |\n'], ...
        k,md_escape(row.id),md_escape(row.status),yes_no(row.solved), ...
        row.time,initial_gap,final_gap,row.iterations,added_cuts, ...
        row.upper_bound,row.lower_bound,row.feasibility_violation);
end
clear cleanup
movefile(temporary,target,'f');
end

function value=get_field(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end

function value=md_escape(value)
value=char(string(value)); value=strrep(value,'|','\|');
value=strrep(value,newline,' ');
end

function value=yes_no(flag)
if flag, value='yes'; else, value='no'; end
end

function files=natural_sort(files)
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
