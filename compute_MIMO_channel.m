function [H, H_I, H_s, E_S] = compute_MIMO_channel(sc_pos, tx_pos, rx_pos, f, a, I_t)
    % compute_MIMO_channel: computes MIMO channels in presence of\ n    % dielectric scatterers using multiple scattering theory.
    %
    % Inputs:
    %   sc_pos  - N×2 scatterer coordinates [x,y]
    %   tx_pos  - M×2 TX coordinates [x,y]
    %   rx_pos  - K×2 RX coordinates [x,y]
    %   f       - frequency (Hz)
    %   a       - N×1 radii of cylinders (m)
    %   I_t     - M×1 TX currents (A)
    %
    % Outputs:
    %   H_I     - incident (LoS) channel, size K×M
    %   H_s     - scattered channel due to cylinders, size K×M
    %   H       - total channel = H_I + H_s, size K×M
    %   E_S     - scattered field at RX terminals, size K×1

    % Physical constants
    mu0   = 4*pi*1e-7;
    eps0  = 8.854e-12;
    omega = 2*pi*f;
    k0    = omega * sqrt(mu0*eps0);

    % Dimensions
    M = size(tx_pos,1);
    K = size(rx_pos,1);
    N = size(sc_pos,1);

    % 1) Incident field matrix (K×M)
    H_I = compute_incident_matrix(rx_pos, tx_pos, k0);

    % 2) Response matrices: (N×M) and (N×K)
    H_t = compute_response_matrix(sc_pos, tx_pos, k0);
    H_r = compute_response_matrix(sc_pos, rx_pos, k0);

    % 3) Multiple scattering: solve (Γ^{-1} - C)^{-1}
    Gamma_inv = compute_Gamma_inv_PEC(k0, a);
    Cmat      = compute_coupling_matrix(sc_pos, k0);
    A         = Gamma_inv - Cmat;
    invA      = A \ eye(N);

    % 4) Scattered channel (K×M)
    H_s = H_r.' * invA * H_t;

    % 5) Total channel and scattered field
    H   = H_I + H_s;
    E_S = -omega*mu0/4 * H_s * I_t;
end
