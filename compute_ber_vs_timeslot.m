function Pe = compute_ber_vs_timeslot(N_obs_total, N, Nc_tx, dt, T_max, Ts_list)
%COMPUTE_BER_VS_TIMESLOT Bit error rate (BER) of the on-off keying S2S
% system as a function of time-slot duration, using the genie-aided
% decision-feedback (DF) detector (Appendix B, Eq. 46-52 of the paper).
%
%   M. Rezaei et al., "Spheroidal Molecular Communication via Diffusion:
%   Signaling Between Homogeneous Cell Aggregates," IEEE Trans. Mol. Biol.
%   Multi-Scale Commun., vol. 10, no. 1, pp. 197-210, Mar. 2024.
%
% METHOD
%   For each candidate time-slot duration Ts:
%     1. Choose the sampling time ts within a slot as the time that
%        maximizes the (normalized) observation probability p_all(t)
%        (Section IV-B: "sampling time is chosen to maximize pobs(t)").
%     2. Determine the number of interfering previous slots, M, that fall
%        within the simulation horizon T_max.
%     3. For every one of the 2^M equiprobable previous-bit patterns
%        (genie-aided: previous bits are assumed perfectly known), compute
%        the Poisson means E0/E1 for the current bit being 0/1 (Eq. 47),
%        the MAP decision threshold (Eq. 49), and the resulting error
%        probability from the Poisson CDFs (Eq. 51-52).
%     4. Average over all previous-bit patterns, weighted equally
%        (Eq. 50), to obtain the overall BER for this Ts.
%
% INPUTS
%   N_obs_total - expected number of observed molecules over time, e.g. the
%                 N_obs_A3 output of COMPUTE_RECEIVER_RESPONSE.m, for a
%                 single transmitted bit "1" (impulse response of the full
%                 S2S system convolved with the release rate)
%   N           - molecules released per cell for bit "1"
%   Nc_tx       - number of cells in the transmitting spheroid
%   dt          - simulation time step corresponding to one sample of
%                 N_obs_total (i.e., 1/w_max from the frequency sweep)
%   T_max       - total simulation time horizon (same units as dt)
%   Ts_list     - vector of candidate time-slot durations to evaluate,
%                 in the same time units as T_max (e.g. seconds)
%
% OUTPUT
%   Pe          - BER for each entry of Ts_list

% Normalized (per-cell, per-molecule) observation probability trajectory.
p_all = N_obs_total / (N * Nc_tx);

% Time index grid, in units of simulation samples (0, 1, 2, ... steps of dt).
tm = 0:1:(T_max/dt);

% Sampling time within a slot: the time index that maximizes p_all
% (skip index 1 / t=0, which is trivially zero for a causal impulse
% response).
[~, nv] = max(p_all(2:end));
ts = tm(nv);

Pe = zeros(1, numel(Ts_list));

for iTs = 1:numel(Ts_list)
    Ts_samples = Ts_list(iTs) / dt;   % slot duration, in simulation samples

    % Number of previous slots that can still interfere within T_max
    % (inter-symbol interference window length).
    M = ceil(T_max / Ts_samples) - 1;

    % Sample times for the current slot (k=0) and the M interfering
    % previous slots (k=1..M), spaced Ts_samples apart (Eq. 17-18).
    sample_times = ts + (0:M) * Ts_samples;
    u = floor(sample_times);
    NN = p_all(u);   % NN(1) = current-slot contribution, NN(2:end) = ISI terms

    Pe_cond = zeros(1, 2^M);

    for w = 0:2^M-1
        % Binary pattern of the M previously-transmitted bits (LSB first).
        % Uses the base-MATLAB BITGET function instead of the Communications
        % Toolbox's dec2binvec, so this runs without extra toolboxes.
        Pre_bits = bitget(w, 1:M);

        isi_term = sum(Pre_bits .* NN(2:end));   % I(ts), Eq. 18

        E0 = Nc_tx * N * isi_term;                       % mean count if bit_0 = 0 (Eq. 47)
        E1 = Nc_tx * N * NN(1) + Nc_tx * N * isi_term;    % mean count if bit_0 = 1 (Eq. 47)

        if E0 == 0
            Thr = 0;
        else
            % MAP decision threshold (Eq. 49).
            Thr = (Nc_tx * N * NN(1)) / log(1 + NN(1) / isi_term);
        end

        % Conditional error probability for this previous-bit pattern
        % (Eq. 51-52), from the Poisson CDF of the two hypotheses.
        Pe_cond(w+1) = 0.5 * (poisscdf(Thr, E1) + 1 - poisscdf(Thr, E0));
    end

    % Average over all 2^M equiprobable previous-bit patterns (Eq. 50).
    Pe(iTs) = (0.5^M) * sum(Pe_cond);
end

end
