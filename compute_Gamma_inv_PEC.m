function Gamma_inv = compute_Gamma_inv_PEC(k0, a)
% compute_Gamma_inv: computes the inverse scattering coefficient matrix for PEC cylinders
% here using PEC boundary condition: Gamma = -J0(k0*a)/H0^{(2)}(k0*a)
% Inputs:
%   k0        - free-space wavenumber
%   epsilon_r - N-by-1 relative permittivities (ignored for PEC)
%   a         cylinder radii
% Output:
%   Gamma_inv - N-by-N diagonal inverse scattering coefficient matrix

N = numel(a);
Gamma_inv = zeros(N);
for i = 1:N
    % PEC scattering for infinite cylinder (m=0 mode)
    gamma_i = -besselj(0, k0*a(i)) / besselh(0, 2, k0*a(i));
    % Invert diagonal element
    Gamma_inv(i,i) = 1 / gamma_i;
end
end
