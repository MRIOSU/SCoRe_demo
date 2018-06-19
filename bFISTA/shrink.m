function w = shrink(s, lam)

T = abs(s);
w = max(T - lam, 0).*s;
T(T == 0) = 1;
w = w./T;