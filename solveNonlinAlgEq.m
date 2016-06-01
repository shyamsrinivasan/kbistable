function solveNonlinAlgEq()

%solve the nonlinear algebraic equation using the SUNDIALS suite of
%nonlinear solvers

%function to be solved i.e f(X) in f(X) = 0 is provideed as an input
%argument to solveNonlinEq()

% Modified/exact Newton
%   msbset = 0   -> modified Newton
%   msbset = 1   -> exact Newton
msbset = 0;

%costraints
%-1 for <= 0
%1 for >= 0
%nconstr = nvar
constraints = [ 0 0 1 -1 1 -1];


options = KINSetOptions('UserData', data, ...
                        'FuncNormTol', ftol, ...
                        'ScaledStepTol', stol, ...
                        'Constraints', constraints, ...
                        'MaxNumSetups', msbset, ...
                        'LinearSolver', 'Dense');

%initialize KINSol
%neq = number of nonlinear equations to solve
KINInit(@f(X),neq,options);

%initial guesses

%scaling of y and f
%yscale
%fscale

%call KINSol
[status, y] = KINSol(y0, strategy, yscale, fscale);
stats = KINGetStats;

%free memory
KINFree;
    