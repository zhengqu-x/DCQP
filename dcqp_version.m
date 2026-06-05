function ver_info = dcqp_version()
% DCQP_VERSION  Return version information for the DC-QP solver
%
% SYNTAX:
%   ver_info = dcqp_version()
%
% OUTPUT:
%   ver_info - Structure with version information:
%              .version     - Version string
%              .date        - Release date
%              .description - Package description

ver_info = struct();
ver_info.version = '1.0.1';
ver_info.date = '2026-06-03';
ver_info.description = 'Doubly Nonnegative based Cutting Plane method for Quadratic Programming';
ver_info.authors = 'Zheng Qu, Defeng Sun, Jintao Xu';
ver_info.license = 'Academic License';
ver_info.url = 'https://github.com/zhengqu-x/DCQP';

if nargout == 0
    % Display version information when called without output
    fprintf('\n=== DCQP Solver ===\n');
    fprintf('Version: %s\n', ver_info.version);
    fprintf('Date: %s\n', ver_info.date);
    fprintf('Description: %s\n', ver_info.description);
    fprintf('License: %s\n', ver_info.license);
    fprintf('URL: %s\n', ver_info.url);
    fprintf('===================\n\n');
    clear ver_info;
end

end
