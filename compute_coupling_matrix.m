function C = compute_coupling_matrix(sc_pos, k0)
% compute_coupling_matrix: N x N mutual coupling via Hankel
N = size(sc_pos,1);
C = zeros(N);
dist = squareform(pdist(sc_pos));
for i=1:N
    for j=1:N
        if i~=j, C(i,j)=besselh(0,2,k0*dist(i,j)); end
    end
end
end