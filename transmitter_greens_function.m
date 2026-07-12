function Conc_tx = transmitter_greens_function(Nc_tx, vc_tx, r_src, theta_src, phi_src, ...
                                                 R_tx, r_obs, theta_obs, phi_obs, ...
                                                 D, w_min, w_max, kd_tx)
%TRANSMITTER_GREENS_FUNCTION Green's function for concentration (GFC) due to
% a single impulsive point-molecule source located INSIDE a homogenized
% porous spheroidal transmitter, evaluated at an arbitrary observation point
% that may lie inside or outside the transmitter.
%
% This implements the boundary-value diffusion problem of Section III-A and
% Appendix A-B (Transmitter Region) of:
%   M. Rezaei et al., "Spheroidal Molecular Communication via Diffusion:
%   Signaling Between Homogeneous Cell Aggregates," IEEE Trans. Mol. Biol.
%   Multi-Scale Commun., vol. 10, no. 1, pp. 197-210, Mar. 2024.
%
% The transmitting spheroid (radius R_tx, effective diffusion coefficient
% D_tx = eps_tx^1.5 * D, Eq. 3) is centered at the coordinate origin. A
% single point source -- representing one releasing cell -- sits at
% (r_src, theta_src, phi_src) inside the spheroid. The concentration is
% solved in the frequency domain (angular frequency omega) as a
% spherical-harmonic series (Eq. 22), with the radial part obtained from the
% linear coefficient system in Eq. (21). The time-domain signal is then
% recovered via an inverse FFT over the sampled frequency grid.
%
% This function is called once per (source-cell, time-grid) combination by
% COMPUTE_RELEASE_RATE.m, which integrates it over the whole transmitter
% volume to obtain the aggregate molecule release rate g(t) (Eq. 10).
%
% INPUTS
%   Nc_tx      - number of cells in the transmitting spheroid
%   vc_tx      - volume of a single cell                              [m^3]
%   r_src, theta_src, phi_src
%              - spherical coordinates of the point source, relative to the
%                spheroid center                                [m, rad, rad]
%   R_tx       - transmitting spheroid (outer) radius                  [m]
%   r_obs, theta_obs, phi_obs
%              - spherical coordinates of the observation point
%   D          - free-medium (extracellular) diffusion coefficient [m^2/s]
%   w_min      - angular-frequency grid step used to build the frequency
%                sweep w_min:w_min:w_max
%   w_max      - upper angular-frequency limit of the sweep
%   kd_tx      - first-order degradation/reaction rate inside the
%                transmitter (K in Eq. 13; set to 0 for a non-reactive
%                transmitter)
%
% OUTPUT
%   Conc_tx    - time-domain concentration response at the observation
%                point, reconstructed on the time grid
%                t = 0 : 1/w_max : (1/domega)*2 + domega, domega = spacing
%                of omega0.

%% Spheroid porosity and effective diffusion coefficient (Eq. 2-3)
eps_tx = 1 - (Nc_tx * vc_tx) / (4*pi/3 * R_tx^3);   % porosity (extracellular fraction)
a_tx   = R_tx;                                       % spheroid outer radius
D_tx   = (eps_tx^1.5) * D;                           % effective diffusion coeff. (tortuosity tau = eps^-0.5, Eq. 3)
alpha_tx = sqrt(D_tx / D);                            % = 1/kappa_tx, used in the boundary-condition rows below

%% Series truncation order
% Number of spherical-harmonic degrees (n = 0..N) retained in series (22).
% N = 4 gives a converged solution for the parameter ranges used in the paper.
N = 4;

%% Spherical Bessel / Hankel functions and their radial derivatives
% jn(.), yn(.): spherical Bessel functions of the first / second kind.
% hn(.): spherical Hankel function of the second kind (the paper's outgoing/
% decaying-wave convention for the unbounded exterior domain).
% The "d" versions are derivatives w.r.t. the function argument, obtained
% from the standard recurrence relation for spherical Bessel functions:
%   f_n'(x) = [n*f_{n-1}(x) - (n+1)*f_{n+1}(x)] / (2n+1)
SphBess   = @(x, n) (pi./(2*x)).^0.5 .* besselj(n+0.5, x);
SphBessn  = @(x, n) (pi./(2*x)).^0.5 .* bessely(n+0.5, x);
SphHank   = @(x, n) (pi./(2*x)).^0.5 .* besselh(n+0.5, 2, x);
SphBessd  = @(x, n) (2*n+1).^-1 .* (n.*SphBess(x, n-1)  - (n+1).*SphBess(x, n+1));
SphBessdn = @(x, n) (2*n+1).^-1 .* (n.*SphBessn(x, n-1) - (n+1).*SphBessn(x, n+1));
SphHankd  = @(x, n) (2*n+1).^-1 .* (n.*SphHank(x, n-1)  - (n+1).*SphHank(x, n+1));

%% Frequency grid and matching time vector for the IFFT reconstruction
omega0 = w_min:w_min:w_max;
domega = omega0(2) - omega0(1);
t = 0:1/max(omega0):(1/domega)*2 + domega;   %#ok<NASGU> % (kept for consistency with callers)

%% Solve the boundary-value problem at every sampled frequency
C = zeros(1, numel(omega0));

for jj = 1:numel(omega0)
    omega = omega0(jj) * pi;

    % Diffusion "wavenumbers" inside (kdp_tx) and outside (kdpo) the
    % spheroid. kdp_tx includes the degradation rate kd_tx (Eq. 13);
    % kdpo is for the reaction-free extracellular medium (Eq. 8).
    kdp_tx = -(kd_tx + 1i*omega) / D_tx;
    kdpo   = -(0     + 1i*omega) / D;

    for l = 0:N
        % Unknown radial-solution coefficients for spherical-harmonic
        % degree l (see Eq. 31-32):
        %   Cn        - amplitude inside the spheroid,  r < r_src
        %   An, Bn    - amplitude inside the spheroid,  r_src < r < a_tx
        %   Dn        - amplitude outside the spheroid, r > a_tx
        %
        % The four rows of CoeffMat enforce, top to bottom:
        %   (1) continuity of concentration at the spheroid boundary a_tx  (Eq. 29)
        %   (2) continuity of diffusive flux at the spheroid boundary a_tx (Eq. 28)
        %   (3) continuity of concentration at the source location r_src   (Eq. 33)
        %   (4) the unit point-source flux jump at r_src                   (Eq. 30)
        CoeffMat = [ ...
            0,                                                          alpha_tx*SphBess(sqrt(kdp_tx)*a_tx, l),                alpha_tx*SphBessn(sqrt(kdp_tx)*a_tx, l),                 -SphHank(sqrt(kdpo)*a_tx, l);
            0,                                                          D_tx*sqrt(kdp_tx)*SphBessd(sqrt(kdp_tx)*a_tx, l),     D_tx*sqrt(kdp_tx)*SphBessdn(sqrt(kdp_tx)*a_tx, l),      -D*sqrt(kdpo)*SphHankd(sqrt(kdpo)*a_tx, l);
            SphBess(sqrt(kdp_tx)*r_src, l),                             -SphBess(sqrt(kdp_tx)*r_src, l),                       -SphBessn(sqrt(kdp_tx)*r_src, l),                        0;
            r_src^2*sqrt(kdp_tx)*SphBessd(sqrt(kdp_tx)*r_src, l),      -r_src^2*sqrt(kdp_tx)*SphBessd(sqrt(kdp_tx)*r_src, l), -r_src^2*sqrt(kdp_tx)*SphBessdn(sqrt(kdp_tx)*r_src, l),   0];

        % Right-hand side: only the flux-jump equation has a non-zero
        % source term, normalized by D_tx (unit-strength point source).
        RHS = [0; 0; 0; 1/D_tx];

        % Solved with the pseudo-inverse for numerical robustness, since
        % CoeffMat can become ill-conditioned for some (l, omega) pairs.
        AS = pinv(CoeffMat) * RHS;

        Cn(l+1) = AS(1); %#ok<AGROW>
        An(l+1) = AS(2); %#ok<AGROW>
        Bn(l+1) = AS(3); %#ok<AGROW>
        Dn(l+1) = AS(4); %#ok<AGROW>
    end

    % Select the radial branch that matches the observation point location.
    ll = 0:N;
    if r_obs <= r_src
        CR = Cn .* SphBess(sqrt(kdp_tx)*r_obs, ll);
    elseif r_obs > a_tx
        CR = Dn .* SphHank(sqrt(kdpo)*r_obs, ll);
    else
        CR = An .* SphBess(sqrt(kdp_tx)*r_obs, ll) + Bn .* SphBessn(sqrt(kdp_tx)*r_obs, ll);
        CR(isnan(CR)) = 0;   % guards against 0*Inf edge cases at r_obs -> 0
    end

    % Assemble the spherical-harmonic sum (Eq. 22-24) for this frequency.
    l1 = 0;
    C(jj) = 0;
    for kt = 1:N+1
        LegP  = legendre(l1, cos(theta_obs));
        LegP0 = legendre(l1, cos(theta_src));

        for m = 0:l1
            Leg  = LegP(m+1);
            Leg0 = LegP0(m+1);

            SphHar  = sqrt((2*l1+1)/(4*pi) * factorial(l1-m)/factorial(l1+m)) * Leg  * cos(m*(phi_obs - phi_src));
            SphHar0 = sqrt((2*l1+1)/(4*pi) * factorial(l1-m)/factorial(l1+m)) * Leg0;

            if m == 0
                kk = 1;   % Eq. 24 normalization: L_0 = 1/(2*pi)
            else
                kk = 2;   % Eq. 24 normalization: L_m = 1/pi, m >= 1
            end

            C(jj) = C(jj) + kk * CR(kt) * SphHar * SphHar0;
        end
        l1 = l1 + 1;
    end
end

%% Reconstruct the time-domain concentration via inverse FFT
% The frequency series C(omega) is Hermitian-extended (conjugate symmetric)
% to build a real-valued time-domain signal.
Conc_tx = real(ifft([0, C*max(omega0), conj(fliplr(C*max(omega0)))]));

end
