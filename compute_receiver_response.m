function [LE, Conc_center_response, Conc_rx_avg, Conc_rx_avg_unscaled, RLS_FUNC, N_obs_A, N_obs_A3, N_obs_A2] = ...
    compute_receiver_response(D, RLS_Rate_total, N, eps_rx, vc_rx, Nc_tx, vc_tx, R_tx, R_rx, r_src, ...
                               center_distance, kd_rx, kd_tx, w_min, w_max)
%COMPUTE_RECEIVER_RESPONSE End-to-end spheroid-to-spheroid (S2S) diffusive
% channel response: convolves the transmitting spheroid's release rate g(t)
% (from COMPUTE_RELEASE_RATE.m) with the receiving spheroid's point-source
% impulse response (RECEIVER_GREENS_FUNCTION.m), volume-averaged over the
% receiver, to obtain the expected number of observed molecules over time
% (Section IV-A, Eq. 14-16 of the paper).
%
%   M. Rezaei et al., "Spheroidal Molecular Communication via Diffusion:
%   Signaling Between Homogeneous Cell Aggregates," IEEE Trans. Mol. Biol.
%   Multi-Scale Commun., vol. 10, no. 1, pp. 197-210, Mar. 2024.
%
% METHOD
%   The transmitter is treated as a point source located at the given
%   center_distance from the receiver (valid when the inter-spheroid
%   distance is much larger than the spheroid sizes themselves, as argued
%   in Section III-A). Its geometric release dynamics are already captured
%   by RLS_Rate_total (the aggregate release rate g(t)).
%   The receiving spheroid is discretized into a 3-D grid of volume
%   elements. For each element, RECEIVER_GREENS_FUNCTION.m gives the
%   impulse response to the point-source transmitter; convolving with
%   RLS_FUNC = N*Nc_tx*RLS_Rate_total gives that element's concentration
%   time series, which is then volume-averaged over the whole receiver
%   (Eq. 14-15) to obtain the expected number of observed molecules (Eq.
%   16).
%
% INPUTS
%   D               - free-medium diffusion coefficient               [m^2/s]
%   RLS_Rate_total  - normalized release rate from COMPUTE_RELEASE_RATE.m
%                     (per one molecule released per cell)
%   N               - molecules released per cell for bit "1"
%   eps_rx          - receiving spheroid porosity (Eq. 2)
%   vc_rx           - volume of a single receiver cell                  [m^3]
%   Nc_tx           - number of cells in the transmitting spheroid
%   vc_tx           - volume of a single transmitter cell (unused directly;
%                     kept for interface consistency with the release-rate
%                     computation that produced RLS_Rate_total)
%   R_tx            - transmitting spheroid radius (unused directly; see
%                     vc_tx note above)                                    [m]
%   R_rx            - receiving spheroid radius                            [m]
%   r_src           - transmitter point-source offset; forced to ~0 inside
%                     this function (point-source approximation for the
%                     already-integrated transmitter geometry)
%   center_distance - center-to-center distance between spheroids           [m]
%   kd_rx           - degradation rate inside the receiver (reaction A->E)
%   kd_tx           - degradation rate inside the transmitter (unused
%                     directly; kept for interface consistency)
%   w_min, w_max    - angular-frequency sweep limits (must match the values
%                     used to generate RLS_Rate_total)
%
% OUTPUTS
%   LE                     - combined signal length (impulse response +
%                            release-rate length - 1), useful for building
%                            a matching time axis for the convolution result
%   Conc_center_response   - concentration response sampled at a single
%                            point near the receiver center (for spot
%                            checks / comparison plots)
%   Conc_rx_avg            - volume-averaged concentration inside the
%                            receiver (Eq. 14-15), spheroid-source case
%   Conc_rx_avg_unscaled    - same quantity before volume normalization by
%                            1/v_rx (diagnostic / intermediate output)
%   RLS_FUNC               - de-normalized release rate used for the
%                            convolutions (N*Nc_tx*RLS_Rate_total)
%   N_obs_A                 - DEPRECATED: kept only for backward interface
%                            compatibility. Flagged as an unreliable
%                            (double-counted) quantity in the original
%                            derivation; not used in the BER analysis.
%   N_obs_A3                - expected number of observed molecules over
%                            time for the full spheroid-transmitter /
%                            spheroid-receiver system (Eq. 16). This is the
%                            primary output used downstream for BER
%                            evaluation.
%   N_obs_A2                - expected number of observed molecules,
%                            rescaled to the "point-source transmitter"
%                            convention (used for the point-source
%                            comparison curves in the paper, e.g. Fig. 8/10)

omega0 = w_min:w_min:w_max;

% Number of cells in the receiving spheroid consistent with the requested
% porosity (inverse of Eq. 2).
Nc_rx = floor((4*pi/3*R_rx^3) * (1 - eps_rx) / vc_rx);

v_rx = 4/3*pi*R_rx^3;   % receiving spheroid volume

% De-normalize the release rate: RLS_Rate_total is the rate per one
% molecule released per cell; scale up to the full transmitter population.
RLS_FUNC = N * Nc_tx .* RLS_Rate_total;

%% Spatial discretization of the receiving spheroid volume
% Volume elements over which the point-source receiver response is
% evaluated and then averaged (Eq. 14-15).
num_r     = 5;
num_theta = 4;
num_phi   = 4;

dr_rx     = R_rx / num_r;
dtheta_rx = pi / num_theta;
dphi_rx   = 2*pi / num_phi;

%% Point-source transmitter approximation
% The transmitter's geometric release dynamics are already folded into
% RLS_FUNC, so here it is treated as a point source essentially at its own
% center (r_src ~ 0), located center_distance away from the receiver.
r_src = 0.001e-6;
phi_src = 0;
theta_src = pi/2;

C_rx_Total  = 0;   % accumulator: volume-weighted convolved response
C_rx_Total2 = 0;   % accumulator: volume-weighted impulse response (no convolution)

for rr = 1:num_r
    r_rx = (rr-1)*dr_rx + dr_rx/2;

    for tt = 1:num_theta
        theta = (tt-1)*dtheta_rx + dtheta_rx/2;

        for ff = 1:num_phi
            phi = (ff-1)*dphi_rx + dphi_rx/2;

            % Impulse response of the receiver at this volume element, due
            % to the point-source transmitter.
            Conc_1_rx = receiver_greens_function(eps_rx, Nc_rx, vc_rx, r_src, theta_src, phi_src, ...
                                                  R_rx, r_rx, theta, phi, ...
                                                  center_distance, D, w_min, w_max, kd_rx);

            % Volume of this (shell, sector) receiver element.
            r1 = (rr-1)*dr_rx;
            r2 = rr*dr_rx;
            V_part = 4/3*pi*(r2^3 - r1^3) * (1/(num_phi*num_theta));

            % Convolve the release rate with this element's impulse
            % response (Eq. 14); the leading 1/(Nc_tx*N) undoes the
            % de-normalization applied to RLS_FUNC above, so this
            % convolution is expressed per one molecule per cell before
            % being re-scaled again for the final outputs below.
            C_RX  = (1/(Nc_tx*N)) .* conv(RLS_FUNC, Conc_1_rx) .* (1/max(omega0));
            C_RX2 = Conc_1_rx;   % raw impulse response, no convolution

            C_rx_Total  = C_rx_Total  + C_RX  .* V_part;
            C_rx_Total2 = C_rx_Total2 + C_RX2 .* V_part;
        end
    end
end

%% Single-point response near the receiver center (diagnostic / plotting)
r_rx_center = 0.0001e-6;
theta_center = pi/2;
phi_center = pi;
Conc_center_impulse = receiver_greens_function(eps_rx, Nc_rx, vc_rx, r_src, theta_src, phi_src, ...
                                                R_rx, r_rx_center, theta_center, phi_center, ...
                                                center_distance, D, w_min, w_max, kd_rx);
Conc_center_response = (1/(Nc_tx*N)) .* conv(RLS_FUNC, Conc_center_impulse) .* (1/max(omega0));

%% Volume-averaged receiver concentration (Eq. 14-15)
Conc_rx_avg          = (1/v_rx) .* C_rx_Total;
Conc_rx_avg_unscaled = (1/v_rx) .* C_rx_Total2;


% Expected number of molecules observed inside the receiver over time for
% the full spheroid transmitter / spheroid receiver system (Eq. 16). This
% is the quantity used downstream for the BER analysis.
P_obs_A2 = v_rx .* Conc_rx_avg_unscaled;
N_obs_A3 = conv(P_obs_A2, RLS_FUNC) .* (1/max(omega0));

% Point-source-transmitter-equivalent scaling (undoes the 1/(Nc_tx*N)
% normalization), used for point-source comparison curves.
N_obs_A2 = P_obs_A2 .* Nc_tx * N;

LE = length(Conc_1_rx) + length(RLS_FUNC) - 1;

end
