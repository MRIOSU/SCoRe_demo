function [uest] = powerIter(A,At,n)

disp('Running power iteration to estimate the Lipschitz constant...');
%Initial guess
q = randn(n);
q = q / norm(q(:));
th = 1e-3; % FISTA code has th = 1e-4, but it's slow
err = inf;
uest = inf;
while err > th
    q = At(A(q));
    %Progress
    unew = norm(q(:));
    err = abs(log(uest / unew));
    %Save the estimate
    uest = unew;
    q = q / norm(q(:));
end
%Use a little more than twice the eigenvalue estimate
% L = 2.05*uest;