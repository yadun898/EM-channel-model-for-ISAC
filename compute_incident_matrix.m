function H = compute_incident_matrix(rx_pos, tx_pos, k0)
    [K,~] = size(rx_pos); [M,~] = size(tx_pos);
    H = zeros(K,M);
    for k = 1:K
        for m = 1:M
            d = norm(rx_pos(k,:) - tx_pos(m,:));
            H(k,m) = besselh(0,2,k0*d);
        end
    end
end