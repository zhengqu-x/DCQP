function watch_structured_summary(result_folder,expected_count,poll_seconds)
%WATCH_STRUCTURED_SUMMARY Compatibility wrapper for WATCH_RESULT_SUMMARY.
if nargin<2, expected_count=[]; end
if nargin<3, poll_seconds=[]; end
watch_result_summary(result_folder,expected_count,poll_seconds);
end
