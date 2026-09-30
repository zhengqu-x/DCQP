function verify_frozen_runtime()
%VERIFY_FROZEN_RUNTIME Verify package integrity and frozen-code provenance.

tests_dir=fileparts(mfilename('fullpath'));
root=fileparts(tests_dir);
packaged=read_manifest(fullfile(tests_dir,'packaged_runtime_sha256.txt'));
frozen=read_manifest(fullfile(tests_dir,'frozen_runtime_sha256.txt'));
approved_reporting_changes={ ...
    'dcqp_check_input.m', ...
    'dcqp_default_params.m', ...
    'dcqp_solve.m', ...
    'dcqp_version.m', ...
    'utils/ExtendedKonnoCut.m', ...
    'utils/ValidatedExtendedKonnoCut.m', ...
    'utils/compute_konno_cut.m', ...
    'utils/compute_lower_bound.m', ...
    'utils/dcqp_konno_step.m', ...
    'utils/dnn_cut_verification.m', ...
    'utils/generate_cut_dnn.m', ...
    'utils/gurobiqp.m', ...
    'utils/lower_bound_dnn.m', ...
    'utils/qpsolver.m', ...
    'utils/refine_failed_cut_point.m'};

for k=1:size(packaged,1)
    expected=packaged{k,1};
    relative=packaged{k,2};
    file=fullfile(root,relative);
    assert(isfile(file),'dcqp:FrozenRuntimeMissing', ...
        'Packaged runtime file is missing: %s',relative);
    actual=sha256_file(file);
    assert(strcmp(actual,expected),'dcqp:FrozenRuntimeChanged', ...
        'Packaged runtime file changed: %s',relative);
end

unchanged=0;
for k=1:size(frozen,1)
    expected=frozen{k,1};
    relative=frozen{k,2};
    if ismember(relative,approved_reporting_changes)
        continue
    end
    assert(strcmp(sha256_file(fullfile(root,relative)),expected), ...
        'dcqp:FrozenRuntimeChanged', ...
        'Unexpected change from frozen runtime: %s',relative);
    unchanged=unchanged+1;
end

fprintf(['Packaged DCQP runtime verified: %d files match; %d remain ', ...
    'byte-identical to the frozen run and %d may contain approved ', ...
    'non-numerical reporting or result-I/O changes.\n'], ...
    size(packaged,1),unchanged,numel(approved_reporting_changes));
end

function entries=read_manifest(file)
lines=splitlines(string(fileread(file)));
lines=lines(strlength(lines)>0);
entries=cell(numel(lines),2);
for k=1:numel(lines)
    parts=split(strtrim(lines(k)));
    entries{k,1}=char(parts(1));
    entries{k,2}=char(parts(end));
end
end

function value=sha256_file(file)
fid=fopen(file,'rb');
assert(fid>=0,'dcqp:FrozenRuntimeIO','Cannot open %s.',file);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
