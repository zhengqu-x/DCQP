function manifest = generate_structured_suite(options)
%GENERATE_STRUCTURED_SUITE Save an inexact StQP and BoxQP test pool.
%   manifest = generate_structured_suite()
%   options.sizes = [12 24 48]; options.seeds = 1:5;
%   options.families = {'stqp_inexact','boxqp_inexact'};
%   options.suite_name = 'pilot_inexact'; options.generator_options = struct();
%   Generator options are nested by family, e.g.
%   options.generator_options.stqp_inexact = struct(...).
%   Writes only below data/structured/generated/<suite_name>. Existing instances are
%   reused only if exactly identical; incompatible files are never overwritten.
if nargin < 1, options = struct(); end
defaults = struct('sizes',[12 24 48],'seeds',1:5, ...
    'families',{{'stqp_inexact','boxqp_inexact'}}, ...
    'suite_name','pilot_inexact','generator_options',struct());
unknown = setdiff(fieldnames(options),fieldnames(defaults));
if ~isempty(unknown), error('dcqp:generationOptions','Unknown option: %s',unknown{1}); end
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(options,fields{k}), options.(fields{k}) = defaults.(fields{k}); end
end
validateattributes(options.sizes,{'numeric'},{'vector','integer','>=',7,'finite','nonempty'});
validateattributes(options.seeds,{'numeric'},{'vector','integer','nonnegative','<=',2^32-1,'nonempty'});
if isstring(options.families), options.families = cellstr(options.families); end
if ~iscellstr(options.families) || isempty(options.families) || ...
        ~all(ismember(options.families,defaults.families))
    error('dcqp:generationFamilies','Unknown or empty benchmark family list.');
end
if numel(unique(options.families)) ~= numel(options.families) || ...
        numel(unique(options.sizes)) ~= numel(options.sizes) || ...
        numel(unique(options.seeds)) ~= numel(options.seeds)
    error('dcqp:generationDuplicates','Families, sizes and seeds must be distinct.');
end
if isempty(regexp(char(options.suite_name),'^[A-Za-z0-9][A-Za-z0-9_-]*$','once'))
    error('dcqp:generationName','suite_name must contain only letters, digits, _ or -.');
end
if ~isstruct(options.generator_options) || ~isscalar(options.generator_options) || ...
        ~all(ismember(fieldnames(options.generator_options),defaults.families))
    error('dcqp:generationOptions','generator_options must be a struct keyed by family.');
end
data_dir = fileparts(fileparts(mfilename('fullpath')));
folder = fullfile(data_dir,'structured','generated',char(options.suite_name));
if ~isfolder(folder), mkdir(folder); end
manifest_file = fullfile(folder,'manifest.mat');
if isfile(manifest_file)
    previous = load(manifest_file,'options');
    if ~isfield(previous,'options') || ~isequaln(previous.options,options)
        error('dcqp:generationOverwrite','Suite configuration changed. Use another suite_name.');
    end
end
manifest = struct('id',{},'family',{},'n',{},'seed',{},'file',{});
for f = 1:numel(options.families)
    family = options.families{f}; generator = str2func(['generate_' family]);
    generator_options = struct();
    if isfield(options.generator_options,family)
        generator_options = options.generator_options.(family);
    end
    for n = options.sizes(:)'
        for seed = options.seeds(:)'
            instance = generator(n,seed,generator_options);
            id = sprintf('%s_n%d_seed%d',family,n,seed);
            file = fullfile(folder,[id '.mat']);
            if isfile(file)
                previous = load(file,'instance');
                if ~isfield(previous,'instance') || ~isequaln(previous.instance,instance)
                    error('dcqp:generationOverwrite','Different instance exists at %s. Use another suite_name.',file);
                end
            else
                save(file,'instance');
            end
            manifest(end+1) = struct('id',id,'family',family,'n',n, ...
                'seed',seed,'file',file); %#ok<AGROW>
        end
    end
end
save(manifest_file,'manifest','options');
writetable(struct2table(manifest,'AsArray',true),fullfile(folder,'manifest.csv'));
fprintf('Saved %d instances in %s\n',numel(manifest),folder);
end
