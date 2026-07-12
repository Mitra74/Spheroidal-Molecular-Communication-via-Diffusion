function RLS_Rate_total = compute_release_rate(D, N, Nc_tx, vc_tx, R_tx, r_src, kd_tx, w_min, w_max)
%COMPUTE_RELEASE_RATE Aggregate molecule release rate g(t) at the boundary
% of a transmitting spheroid, given an instantaneous release of N molecules
% per cell at t = 0 (Section III-A, Eq. 9-10 of the paper).
%
%   M. Rezaei et al., "Spheroidal Molecular Communication via Diffusion:
%   Signaling Between Homogeneous Cell Aggregates," IEEE Trans. Mol. Biol.
%   Multi-Scale Commun., vol. 10, no. 1, pp. 197-210, Mar. 2024.
%
% METHOD
%   1. The transmitting spheroid is discretized into a 3-D grid of volume
%      elements (num_r x num_theta x num_phi). Each element is treated as a
%      point source and its contribution to the concentration OUTSIDE the
%      spheroid is obtained from TRANSMITTER_GREENS_FUNCTION.m.
%   2. Per-cell contributions are volume-weighted and summed to give the
%      total concentration just outside the spheroid, at a set of
%      concentric spherical shells (radial "bins") extending outward from
%      the spheroid surface.
%   3. Those shell concentrations are volume-integrated to obtain the
%      total number of molecules that have left the spheroid, Nout(t)
%      (unnormalized here; normalization by Nc_tx*N happens in the caller).
%   4. The release rate is the time derivative of Nout(t) (Eq. 10),
%      approximated here by a first-order finite difference.
%
% NOTE: this is the computationally expensive step of the pipeline (nested
% loops over spatial grid x frequency sweep x Bessel-function evaluations).
% It is recommended to run it once per (Nc_tx, kd_tx) combination and cache
% the result (see MAIN_S2S_BER_ANALYSIS.m, which loads/saves a .mat file).
%
% INPUTS
%   D        - free-medium diffusion coefficient                   [m^2/s]
%   N        - molecules released per cell for bit "1"
%   Nc_tx    - number of cells in the transmitting spheroid
%   vc_tx    - volume of a single cell                                [m^3]
%   R_tx     - transmitting spheroid radius                            [m]
%   r_src    - placeholder source radial coordinate, overwritten inside the
%              spatial-discretization loop below (kept for signature
%              compatibility with earlier versions of this function)
%   kd_tx    - first-order degradation rate inside the transmitter (0 for a
%              non-reactive transmitter)
%   w_min    - angular-frequency grid step
%   w_max    - upper angular-frequency limit
%
% OUTPUT
%   RLS_Rate_total - discrete-time samples of the release rate g(t) (per
%                    unit "one molecule per cell"), on the time grid implied
%                    by w_min, w_max (spacing = 1/w_max).

omega0 = w_min:w_min:w_max;
domega = omega0(2) - omega0(1);
t = 0:1/max(omega0):(1/domega)*2 + domega;

v_tx = 4/3*pi*R_tx^3;   % transmitting spheroid volume

%% Source point (varies over the discretization loop below)
phi_src = 0;
theta_src = pi/2;

%% Observation point (fixed direction; only the radial distance varies)
phi_obs = 0;
theta_obs = pi/2;

%% Spatial discretization of the transmitting spheroid volume
% The spheroid is divided into num_r radial shells, each split into
% num_theta x num_phi angular sectors. Each (shell, sector) cell is
% represented by a point source at its centroid.
num_r     = 6;
num_theta = 4;
num_phi   = 4;

dr_tx     = R_tx / num_r;
dtheta_tx = pi / num_theta;
dphi_tx   = 2*pi / num_phi;

%% Radial "shells" outside the spheroid used to accumulate Nout(t)
% 20 shells of thickness R_tx starting at the spheroid surface, i.e. the
% concentration is sampled well beyond the transmitter to capture the
% outward-diffusing population.
Conc_TOTAL_OUT = 0;
v_total_OUT = 0;

for robs = 1:20
    dr_obs = 1*R_tx;
    r_obs  = R_tx + (robs-1)*dr_obs + dr_obs/2;   % shell centroid (observation radius)
    r1_out = R_tx + (robs-1)*dr_obs;              % shell inner radius
    r2_out = R_tx + robs*dr_obs;                  % shell outer radius

    Conc_Total = 0;

    for rr = 1:num_r
        r_src = (rr-1)*dr_tx + dr_tx/2;

        for tt = 1:num_theta
            theta_src = (tt-1)*dtheta_tx + dtheta_tx/2;

            for ff = 1:num_phi
                phi_src = (ff-1)*dphi_tx + dphi_tx/2;

                % Concentration outside the spheroid due to this single
                % source point, evaluated at the current shell radius.
                Conc_1_cell = transmitter_greens_function(Nc_tx, vc_tx, r_src, theta_src, phi_src, ...
                                                           R_tx, r_obs, theta_obs, phi_obs, ...
                                                           D, w_min, w_max, kd_tx);
                Conc_1_cell = movmean(Conc_1_cell, 1);   % (no-op smoothing; kept for reproducibility)

                % Volume of this (shell, sector) source element.
                r1 = (rr-1)*dr_tx;
                r2 = rr*dr_tx;
                V_part = 4/3*pi*(r2^3 - r1^3) * (1/(num_phi*num_theta));

                Conc_Total = Conc_Total + Conc_1_cell .* V_part;
            end
        end
    end

    Conc_Total = (1/v_tx) .* Conc_Total;   % normalize by spheroid volume
    Conc_Total(1:30) = max(Conc_Total(1:30), 0);   % clip small negative numerical noise near t = 0

    % Integrate this shell's concentration over its own volume and
    % accumulate into the total "outside" population Nout(t).
    V_part_OUT = 4/3*pi*(r2_out^3 - r1_out^3);
    v_total_OUT = v_total_OUT + V_part_OUT;
    Conc_TOTAL_OUT = Conc_TOTAL_OUT + Conc_Total .* V_part_OUT;
end

%% Differentiate Nout(t) to obtain the release rate g(t), Eq. (10)
% The IFFT-based reconstruction is only reliable up to the point where the
% (numerically noisy) tail starts to turn non-physical (Nout should be
% monotonically non-decreasing); truncate there before differentiating.
idx = find(diff(Conc_TOTAL_OUT) < 0, 1);
C_cut = Conc_TOTAL_OUT(1:idx);

dt = t(2) - t(1);
RLS_Rate_total = diff(C_cut) ./ dt;

end
