function H = compute_response_matrix(sc_pos, ent_pos, k0)

N = size(sc_pos,1);
E = size(ent_pos,1);
H = zeros(N,E);
for i=1:N
    for e=1:E
        d = norm(sc_pos(i,:) - ent_pos(e,:));
        H(i,e) = besselh(0,2,k0*d);
    end
end
end