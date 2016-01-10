% Work out the number of bits used for BLAS/LAPACK indexes. Using the right code
% and/or compiler options, depending on blas_64bit_ints, is really important!
[dummy, max_array_length] = computer;
blas_64bit_ints = max_array_length > 2^31;

% Linking options for BLAS/LAPACK
opts = {};
if exist('octave_config_info')
    % Octave
    link_opts = {'-llapack', '-lblas'};
else
    % Matlab
    link_opts = {'-lmwlapack', '-lmwblas'};
    if blas_64bit_ints
        opts = {'-largeArrayDims'};
    end
end

% While Octave and recent Matlab can compile everything on one step, older
% versions of Matlab refused to be passed .f and .c files in one mex call.

% So, first compile the appropriate Fortran code
% An alternative to selecting different .f files is to set compiler options.
% See ../make_ilp64.sh for a discussion, and a commented out alternative at the
% bottom of this file.
basenames = {'dpofrt', 'dpo2ft'};
if blas_64bit_ints
    suffix = '_ilp64';
else
    suffix = '';
end
files = cellfun(@(x) ['..', filesep, x, suffix, '.f'], ...
        basenames, 'UniformOutput', 0);
mex('-c', files{:});

% Then link everything together
if ispc
    obj = '.obj';
else
    obj = '.o';
end
obj_files = cellfun(@(x) [x, suffix, obj], basenames, 'UniformOutput', 0);
if blas_64bit_ints
    opts = {opts{:}, '-DBLAS64INT'};
end
mex(opts{:}, 'chol_rev.c', obj_files{:}, link_opts{:})

% Clean up
cellfun(@(x) delete(x), obj_files);
if exist('octave_config_info')
    unlink('chol_rev.o'); % only Octave leaves this file behind
end


% Under linux I used to use the procedure below, which can give faster binaries,
% but relies on knowing the compiler. The correct solution is probably for users
% to add desirable optimizations like "-march=native" to their local mex
% configuration.
%
%if blas_64bit_ints
%    system('gfortran -fdefault-integer-8 -march=native -Wall -Wextra -O3 -fPIC -c ../dpo2ft.f');
%    system('gfortran -fdefault-integer-8 -march=native -Wall -Wextra -O3 -fPIC -c ../dpofrt.f');
%else
%    system('gfortran -march=native -Wall -Wextra -O3 -fPIC -c ../dpo2ft.f');
%    system('gfortran -march=native -Wall -Wextra -O3 -fPIC -c ../dpofrt.f');
%end
%mex('chol_rev.c', 'dpofrt.o', 'dpo2ft.o', opts{:}, link_opts{:})

