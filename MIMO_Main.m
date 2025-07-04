%% Fixed TX/RX Array and Single-Target Scattering Channel Simulation
clear; clc;

%% 1) Load configuration parameters into Cfg struct
Cfg.f           = 20e9;                                      % Simulation frequency (Hz)
Cfg.c           = physconst('LightSpeed');                  % Speed of light (m/s)
Cfg.lambda0     = Cfg.c / Cfg.f;                            % Wavelength (m)
Cfg.scalingparam= 0.1;                                       % Scatterer size factor


%% 2) Define fixed TX and RX arrays using helper function
% Example parameters for ULAs
d_center_rx = [10,0];         % Center of RX array
tilt_rx     = pi/2;            % Rotation angle (rad)
N_rx        = 8;            % Number of RX elements
d_rx        = 0.5*Cfg.lambda0; % Spacing

d_center_tx = [0,0];
tilt_tx     = pi/2;         
N_tx        = 8;
d_tx        = 0.5*Cfg.lambda0;

% Generate arrays
rx_array = generate_linear_array(d_center_rx, tilt_rx, N_rx, d_rx);
tx_array = generate_linear_array(d_center_tx, tilt_tx, N_tx, d_tx);

I_t = ones(N_tx,1);                                % Transmit current element

%% 3) Define single-target scatterer positions using equilateral triangle generator
% Scatterers along edges of an equilateral triangle centered at target center
sideLength = 1.0 * Cfg.lambda0;                           % Side length
orient = pi/6;                                           % Triangle orientation (rad)
sc_center = [0.5*Cfg.lambda0, 0.5*Cfg.lambda0];        % Triangle center
sc_pos = generate_equilateral_triangle_scatterers(sideLength, orient, sc_center, Cfg.lambda0, Cfg.scalingparam);

% Display scatterer positions
disp(sc_pos);


%% 4) Compute MIMO channel matrix H and scattered field E_S
[H, ~, ~, E_S] = compute_MIMO_channel(sc_pos, tx_array, rx_array, ...
                    Cfg.f, repmat(Cfg.scalingparam, size(sc_pos,1),1), ...
                    I_t);

%% 5) Save and display results
save('fixed_mimo_channel.mat', 'rx_array','tx_array','sc_pos','H','E_S');
disp('Fixed-array MIMO channel and scattered field computed and saved.');

%% Helper function: generate_linear_array
function pos = generate_linear_array(center, tilt, N, spacing)
    coords = ((0:N-1)' - (N-1)/2) * spacing;
    positions = [coords, zeros(N,1)];
    R = [cos(tilt), -sin(tilt); sin(tilt), cos(tilt)];
    pos = (R * positions')' + center;
end

%% Helper function: generate_equilateral_triangle_scatterers
function sc_pos = generate_equilateral_triangle_scatterers(sideLength, orient, center, wavelength, spacing_factor)
    % Generate scatterer positions along edges of an equilateral triangle
    % sideLength: length of each edge
    % orient: rotation angle (rad)
    % center: [x,y] translation
    % wavelength: for computing spacing
    % spacing_factor: multiplier for spacing
    if nargin < 5
        spacing_factor = 0.1;
    end
    spacing = spacing_factor * wavelength;

    % Define equilateral triangle vertices
    h = sideLength * sqrt(3)/2;
    verts = [0, 0; sideLength, 0; sideLength/2, h];

    % Center at origin
    verts = verts - mean(verts,1);

    % Sample points along edges
    sc_pos = [];
    for k = 1:3
        p1 = verts(k,:);
        p2 = verts(mod(k,3)+1,:);
        d = norm(p2-p1);
        npts = max(ceil(d/spacing)+1,2);
        t = linspace(0,1,npts)';
        pts = p1 + (p2-p1).*t;
        if k>1, pts(1,:)=[]; end
        sc_pos = [sc_pos; pts];
    end

    % Rotate and translate
    Rmat = [cos(orient), -sin(orient); sin(orient), cos(orient)];
    sc_pos = (Rmat*sc_pos')' + center;
    sc_pos = unique(sc_pos,'rows','stable');
end
