function manifest = prepare_structured_pilot()
%PREPARE_STRUCTURED_PILOT Prepare the 30 active cases for equality-free DCQP.
%   StQPs eliminate their last variable; BoxQPs retain their original model.
%   Originals, affine maps, objective offsets and validation are saved with
%   each prepared instance. The source pilot files are never modified.
data_dir = fileparts(fileparts(mfilename('fullpath')));
source = generate_structured_suite();
[~,order] = sortrows([[source.n]',[source.seed]'],[1 2]);
source = source(order);
folder = fullfile(data_dir,'structured','generated','pilot_inexact_reduced');
if ~isfolder(folder), mkdir(folder); end
manifest = source;
for k = 1:numel(source)
    loaded = load(source(k).file,'instance'); original = loaded.instance;
    n = size(original.Q,1);
    if isempty(original.Aeq)
        instance = original;
        E = eye(n); h = zeros(n,1); m = n;
        objective_constant = 0; expansion_error = 0;
        representation = 'original_box';
        id = source(k).id;
    else
        assert(isequal(original.Aeq,ones(1,n)) && isequal(original.beq,1));
        assert(isequal(original.A,[-eye(n);eye(n)]) && ...
            isequal(original.b,[zeros(n,1);ones(n,1)]));
        m = n-1; E = [eye(m);-ones(1,m)]; h = [zeros(m,1);1];
        Q = E'*original.Q*E; Q = (Q+Q')/2;
        c = full(E'*(2*original.Q*h+original.c));
        objective_constant = h'*original.Q*h+original.c'*h;
        instance = struct('Q',Q,'c',c,'A',[-eye(m);ones(1,m)], ...
            'b',[zeros(m,1);1],'Aeq',zeros(0,m),'beq',zeros(0,1), ...
            'LB',zeros(m,1),'UB',ones(m,1),'metadata',original.metadata);
        instance.metadata.known_optimum = original.metadata.known_optimum-objective_constant;
        instance.metadata.known_solution = original.metadata.known_solution(1:m);
        if isfield(instance.metadata,'known_nonglobal_solution')
            point = instance.metadata.known_nonglobal_solution;
            if ~isempty(point), instance.metadata.known_nonglobal_solution = point(1:m); end
            instance.metadata.known_nonglobal_value = ...
                instance.metadata.known_nonglobal_value-objective_constant;
        end
        instance.metadata.objective_convention = ...
            'y''*Q*y+c''*y; add objective_constant for the original objective';
        representation = 'last_simplex_variable_eliminated';
        id = [source(k).id '_reduced'];

        % Verify all coefficients through vertices and edge midpoints, plus
        % the centroid and a known minimizer. All points remain feasible.
        vertices = [zeros(m,1),eye(m)];
        Y = [vertices,ones(m,1)/n,instance.metadata.known_solution];
        for i = 1:n
            for j = i+1:n
                Y(:,end+1) = (vertices(:,i)+vertices(:,j))/2; %#ok<AGROW>
            end
        end
        X = E*Y+h;
        full_values = sum(X.*(original.Q*X),1)+original.c'*X;
        reduced_values = sum(Y.*(Q*Y),1)+c'*Y+objective_constant;
        expansion_error = max(abs(full_values-reduced_values));
        assert(expansion_error<=1e-12*max(1,norm(original.Q,Inf)));
        assert(max(max(original.A*X-original.b))<=1e-12);
        assert(max(abs(original.Aeq*X-original.beq))<=1e-12);
    end
    instance.metadata.source_id = source(k).id;
    instance.metadata.original_n = n;
    instance.metadata.objective_constant = objective_constant;
    instance.metadata.representation = representation;
    transform = struct('E',E,'h',h,'objective_constant',objective_constant, ...
        'source_id',source(k).id,'source_file',source(k).file,'original_n',n);
    file = fullfile(folder,[id '.mat']);
    if isfile(file)
        previous = load(file,'instance','original','transform');
        assert(isequaln(previous.instance,instance) && ...
            isequaln(previous.original,original) && isequaln(previous.transform,transform), ...
            'Existing prepared instance differs; preserve it and choose another folder.');
    else
        save(file,'instance','original','transform','expansion_error');
    end
    manifest(k) = struct('id',id,'family',source(k).family,'n',m, ...
        'seed',source(k).seed,'file',file);
    fprintf('%2d/%d %-38s objective expansion error %.3e\n', ...
        k,numel(source),id,expansion_error);
end
save(fullfile(folder,'manifest.mat'),'manifest','source');
writetable(struct2table(manifest),fullfile(folder,'manifest.csv'));
fprintf('Prepared %d cases in %s\n',numel(manifest),folder);
end
