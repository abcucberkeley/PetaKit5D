function [] = atCall(counterFn, n, fcn)
% counts its calls in the file counterFn and at the n-th runs fcn (n a
% vector and fcn a cell array: at the n(i)-th, fcn{i}; in order if several
% have the same n), e.g., as the call of a pktest.FunctionSpyFixture of
% pause, to make the files a job would make appear during a wait

k = 1;
if isfile(counterFn)
    k = str2double(fileread(counterFn)) + 1;
end
fid = fopen(counterFn, 'w');
fprintf(fid, '%d', k);
fclose(fid);
if ~iscell(fcn)
    fcn = {fcn};
end
for i = reshape(find(k == n), 1, [])
    fcn{i}();
end

end
