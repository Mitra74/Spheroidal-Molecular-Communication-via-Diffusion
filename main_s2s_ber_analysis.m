%% MAIN_S2S_BER_ANALYSIS
% End-to-end spheroid-to-spheroid (S2S) molecular communication example.
%
% Reproduces (a simplified, single-parameter-set version of) the type of
% results shown in Figs. 6-9 of:
%   M. Rezaei et al., "Spheroidal Molecular Communication via Diffusion:
%   Signaling Between Homogeneous Cell Aggregates," IEEE Trans. Mol. Biol.
%   Multi-Scale Commun., vol. 10, no. 1, pp. 197-210, Mar. 2024.
%
% PIPELINE
%   1. Compute (or load a cached) release rate g(t) for the transmitting
%      spheroid -> src/compute_release_rate.m
%   2. For each receiving-spheroid porosity (plus a "transparent receiver"
%      / free-space reference case), compute the expected number of
%      molecules observed inside the receiver over time
%      -> src/compute_receiver_response.m
%   3. For each of those cases, evaluate the bit error rate (BER) of an
%      on-off keying system across a range of time-slot durations
%      -> src/compute_ber_vs_timeslot.m
%

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));

%% ---- Physical / geometric parameters (Table I of the paper) ----
D  = 1e-9;          % free-medium diffusion coefficient                [m^2/s]
N  = 1;             % molecules released per cell for bit "1"
distance = 1000e-6; % center-to-center distance between spheroids         [m]

R_tx = 275e-6;      % transmitting spheroid radius                        [m]
R_rx = 275e-6;      % receiving spheroid radius                           [m]
vc_tx = 3.14e-15;   % single transmitter-cell volume (HepaRG cell)      [m^3]
vc_rx = 3.14e-15;   % single receiver-cell volume                       [m^3]
Nc_tx = 24000;      % number of cells in the transmitting spheroid (fixed)

kd_tx = 0;          % transmitter degradation rate (non-reactive transmitter)
kd_rx = 0.01;       % receiver degradation rate (reaction A -> E, Eq. 1)

r_src = 0.001e-6;   % point-source placeholder offset (see compute_release_rate.m)

% Receiving-spheroid porosities to compare (Eq. 2). eps = 0.1349 corresponds to Nc = 24000 HepaRG cells.
EPS_RX_LIST = [0.1349, 0.2791, 0.4233, 0.5675];

%% ---- Frequency sweep / time grid used throughout the pipeline ----
w_min = 1e-6;   % angular-frequency grid step
w_max = 0.1;    % upper angular-frequency limit
omega0 = w_min:w_min:w_max;
domega = omega0(2) - omega0(1);
t = 0:1/max(omega0):(1/domega)*2 + domega;

dt = 1/max(omega0);        % simulation time step
T_max = 12000/dt;          % simulation horizon for the BER analysis

v_rx = 4/3*pi*R_rx^3;

%% ---- Step 1: transmitting-spheroid release rate g(t) ----
% This is the most expensive step (nested spatial x frequency loops), so
% the result is cached to disk and reused across runs.
release_rate_cache = fullfile(fileparts(mfilename('fullpath')), 'data', ...
                               sprintf('release_rate_Nc%d.mat', Nc_tx));

if exist(release_rate_cache, 'file')
    fprintf('Loading cached release rate from %s\n', release_rate_cache);
    loaded = load(release_rate_cache, 'RLS_Rate_total');
    RLS_Rate_total = loaded.RLS_Rate_total;
else
    fprintf('Computing release rate g(t) for Nc_tx = %d (this can take a while)...\n', Nc_tx);
    RLS_Rate_total = compute_release_rate(D, N, Nc_tx, vc_tx, R_tx, r_src, kd_tx, w_min, w_max);
    save(release_rate_cache, 'RLS_Rate_total');
end

RLS_FUNC = N .* Nc_tx .* RLS_Rate_total;

%% ---- Reference case: free-space ("transparent receiver") point source ----
% Classical 3-D diffusion Green's function for an instantaneous point
% source with first-order degradation (no porous receiver structure). Used
% as the non-porous baseline comparison in Figs. 7-10 of the paper.
CPT = 1 ./ ((4*pi*D*t).^1.5) .* exp(-kd_rx*t - distance^2 ./ (4*D*t));
CPT(1) = 0;   % avoid the t=0 singularity

%% ---- Step 2 & 3: concentration + BER for each receiver configuration ----
n_cases = numel(EPS_RX_LIST) + 1;   % + 1 for the transparent-receiver case
Ts_list = 600:50:1200;              % candidate time-slot durations [s]

N_obs_all = cell(1, n_cases);
Pe_all    = cell(1, n_cases);
legend_labels = cell(1, n_cases);

for case_idx = 1:n_cases
    if case_idx <= numel(EPS_RX_LIST)
        %% Porous spheroidal receiver
        eps_rx = EPS_RX_LIST(case_idx);
        fprintf('Case %d/%d: spheroidal receiver, eps_rx = %.4f\n', case_idx, n_cases, eps_rx);

        [~, ~, ~, ~, ~, ~, N_obs_A3, ~] = compute_receiver_response( ...
            D, RLS_Rate_total, N, eps_rx, vc_rx, Nc_tx, vc_tx, R_tx, R_rx, r_src, ...
            distance, kd_rx, kd_tx, w_min, w_max);

        legend_labels{case_idx} = sprintf('Spheroidal receiver, \\epsilon_{rx} = %.4f', eps_rx);
    else
        %% Transparent (non-porous, free-space) receiver baseline
        fprintf('Case %d/%d: transparent receiver (free-space baseline)\n', case_idx, n_cases);

        P_obs_A = v_rx .* CPT;
        N_obs_A3 = conv(P_obs_A, RLS_FUNC) .* (1/max(omega0));

        legend_labels{case_idx} = 'Transparent receiver';
    end

    N_obs_all{case_idx} = N_obs_A3;
    Pe_all{case_idx} = compute_ber_vs_timeslot(N_obs_A3, N, Nc_tx, dt, T_max, Ts_list);
end

%% ---- Plot: expected number of observed molecules over time ----
figure('Name', 'Expected observed molecules');
hold on;
for case_idx = 1:n_cases
    plot(t(1:1000), N_obs_all{case_idx}(1:1000));
end
xlabel('Time [s]');
ylabel('Expected number of observed molecules');
title('Receiver response for different porosities (Eq. 14-16)');
legend(legend_labels, 'Location', 'best');

%% ---- Plot: BER vs. time-slot duration ----
figure('Name', 'BER vs time slot');
hold on;
for case_idx = 1:n_cases
    semilogy(Ts_list, Pe_all{case_idx});
end
set(gca, 'YScale', 'log');
xlabel('Time slot duration T_s [s]');
ylabel('Bit error rate');
title('BER vs. time-slot duration (Appendix B, Eq. 46-52)');
legend(legend_labels, 'Location', 'best');
