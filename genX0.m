function x0 = genX0(n,R,T,NT,D,SS,lb,ub)

x0 = nan(n,numel(lb));
for i = 1:n
    fin = 0;
    while ~fin
        x =  rand(numel(lb),1) .* (ub-lb) + lb;
        ll = dpLL(R, T, NT, D, SS, x);
        if ~isinf(ll)
            fin = 1;
        end
    end
    x0(i,:) = x;
end