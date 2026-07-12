function Conc_rx = receiver_greens_function(eps_rx, Nc_rx, vc_rx, r_src, theta_src, phi_src, ...
                                              R_rx, r_obs, theta_obs, phi_obs, ...
                                              center_distance, D, w_min, w_max, kd_rx)
%RECEIVER_GREENS_FUNCTION Green's function for concentration (GFC) due to a
% single impulsive point-molecule source located OUTSIDE a homogenized
% porous spheroidal receiver, evaluated at an arbitrary observation point
% that may lie inside or outside the receiver.
%
% This implements the boundary-value diffusion problem of Section III-B and
% Appendix A-C (Receiver Region) of:
%   M. Rezaei et al., "Spheroidal Molecular Communication via Diffusion:
%   Signaling Between Homogeneous Cell Aggregates," IEEE Trans. Mol. Biol.
%   Multi-Scale Commun., vol. 10, no. 1, pp. 197-210, Mar. 2024.
%
% The receiving spheroid (radius R_rx, effective diffusion coefficient
% D_rx = eps_rx^1.5 * D, Eq. 3) is centered at the coordinate origin. A
% single point source sits OUTSIDE the spheroid; its radial coordinate is
% computed here as the distance between the two spheroid centers minus the
% source's offset within the transmitting spheroid (r_src), which is how
% this function is invoked from COMPUTE_RECEIVER_RESPONSE.m for a
% transmitter located a fixed "center_distance" away from the receiver.
%
% Unlike the transmitter case (source inside), here the interior of the
% spheroid is a single radial branch (Gn/Cn from r = 0 to r = R_rx), while
% the exterior is split into a near branch (An, Bn, between the spheroid
% boundary and the source) and a far branch (Dn, beyond the source),
% matching Eq. (43)-(44) and the linear system in Eq. (36).
%
% INPUTS
%   eps_rx      - porosity of the receiving spheroid (Eq. 2)
%   Nc_rx       - number of cells in the receiving spheroid (accepted for
%                 interface consistency with the transmitter-side function;
%                 not used directly since eps_rx already encodes porosity)
%   vc_rx       - volume of a single receiver cell [m^3] (see note above;
%                 not used directly)
%   r_src, theta_src, phi_src
%               - spherical coordinates of the point source's offset within
%                 the transmitting spheroid [m, rad, rad]
%   R_rx        - receiving spheroid (outer) radius [m]
%   r_obs, theta_obs, phi_obs
%               - spherical coordinates of the observation point inside/
%                 outside the receiver
%   center_distance
%               - center-to-center distance between transmitting and
%                 receiving spheroids [m]
%   D           - free-medium (extracellular) diffusion coefficient [m^2/s]
%   w_min       - angular-frequency grid step (sweep w_min:w_min:w_max)
%   w_max       - upper angular-frequency limit of the sweep
%   kd_rx       - first-order degradation/reaction rate inside the
%                 receiver (K in Eq. 13; the reaction A -> E of Eq. 1)
%
% OUTPUT
%   Conc_rx     - time-domain concentration response at the observation
%                 point.

%% Effective radial position of the point source relative to the receiver
% The source lies at radial offset r_src inside the transmitter, which is
% itself centered center_distance away from the receiver (co-linear
% arrangement, consistent with theta_src = pi/2, phi_src = 0 in the caller).
r_src = abs(center_distance - r_src);

a_rx = R_rx;                       % receiving spheroid outer radius
D_rx = (eps_rx^1.5) * D;           % effective diffusion coeff. inside the receiver (Eq. 3)
alpha_rx = sqrt(D_rx / D);         % = 1/kappa_rx, used in the boundary-condition rows below

%% Series truncation order (spherical-harmonic degrees n = 0..N, Eq. 22)
N = 4;

%% Spherical Bessel / Hankel functions and their radial derivatives
% (see TRANSMITTER_GREENS_FUNCTION.m for the definitions/recurrence used)
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

    % Diffusion "wavenumbers" inside (kdp_rx, includes degradation kd_rx,
    % Eq. 13) and outside (kdpo, reaction-free medium, Eq. 12).
    kdp_rx = -(kd_rx + 1i*omega) / D_rx;
    kdpo   = -(0     + 1i*omega) / D;

    for l = 0:N
        % Unknown radial-solution coefficients for spherical-harmonic
        % degree l (see Eq. 43-44):
        %   Gn (Cn)   - amplitude inside the spheroid,      r < a_rx
        %   An, Bn    - amplitude outside the spheroid,     a_rx < r < r_src
        %   Dn        - amplitude outside the spheroid,     r > r_src
        %
        % The four rows of CoeffMat enforce, top to bottom:
        %   (1) continuity of concentration at the spheroid boundary a_rx  (Eq. 41)
        %   (2) continuity of diffusive flux at the spheroid boundary a_rx (Eq. 40)
        %   (3) continuity of concentration at the source location r_src   (Eq. 45)
        %   (4) the unit point-source flux jump at r_src                   (Eq. 42)
        CoeffMat = [ ...
            alpha_rx*SphBess(sqrt(kdp_rx)*a_rx, l),                     -SphBess(sqrt(kdpo)*a_rx, l),                      -SphBessn(sqrt(kdpo)*a_rx, l),                        0;
            D_rx*sqrt(kdp_rx)*SphBessd(sqrt(kdp_rx)*a_rx, l),           -D*sqrt(kdpo)*SphBessd(sqrt(kdpo)*a_rx, l),        -D*sqrt(kdpo)*SphBessdn(sqrt(kdpo)*a_rx, l),          0;
            0,                                                          SphBess(sqrt(kdpo)*r_src, l),                      SphBessn(sqrt(kdpo)*r_src, l),                       -SphHank(sqrt(kdpo)*r_src, l);
            0,                                                          r_src^2*sqrt(kdpo)*SphBessd(sqrt(kdpo)*r_src, l),  r_src^2*sqrt(kdpo)*SphBessdn(sqrt(kdpo)*r_src, l),   -r_src^2*sqrt(kdpo)*SphHankd(sqrt(kdpo)*r_src, l)];

        % Right-hand side: only the flux-jump equation has a non-zero
        % source term, normalized by D (unit-strength point source in the
        % free-diffusion exterior medium).
        RHS = [0; 0; 0; 1/D];

        % Solved with the pseudo-inverse for numerical robustness.
        AS = pinv(CoeffMat) * RHS;

        Cn(l+1) = AS(1); %#ok<AGROW>
        An(l+1) = AS(2); %#ok<AGROW>
        Bn(l+1) = AS(3); %#ok<AGROW>
        Dn(l+1) = AS(4); %#ok<AGROW>
    end

    % Select the radial branch that matches the observation point location.
    ll = 0:N;
    if r_obs <= a_rx
        CR = Cn .* SphBess(sqrt(kdp_rx)*r_obs, ll);
    elseif r_obs > r_src
        CR = Dn .* SphHank(sqrt(kdpo)*r_obs, ll);
    else
        CR = An .* SphBess(sqrt(kdpo)*r_obs, ll) + Bn .* SphBessn(sqrt(kdpo)*r_obs, ll);
        CR(isnan(CR)) = 0;   % guards against 0*Inf edge cases
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
Conc_rx = real(ifft([0, C*max(omega0), conj(fliplr(C*max(omega0)))]));

end
