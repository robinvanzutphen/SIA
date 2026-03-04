%% ========================================================================
%  DEMO 6 — Overlapping Gaussian source + Gaussian detector (same center)
%           p_eff(rho) for non-uniform illumination and detection
%
%
%  We compare 5 curves:
%   1) Sampled p_eff(rho)  (importance sampling from the weighted 2D PDFs)
%   2) Baseline p_eff(rho) (real-space distance-averaging quadrature; stable,
%                           directly matches the manuscript delta-distance definition)
%   3) Hankel p_eff(rho)   (Fourier–Bessel / J0 integral; cross-check)
%   4) Naive: p_unif(rho) * exp(-2 rho^2 / w^2)  (renormalized)
%   5) Naive: p_unif(rho) * exp(-4 rho^2 / w^2)  (renormalized)
%
%  ------------------------------------------------------------------------
%  Manuscript definition (what we approximate numerically)
%  ------------------------------------------------------------------------
%  For concentric / overlapping apertures, define normalized 2D densities
%  q_s(rs) and q_d(rd) over the surface (units 1/area), such that:
%     ∬_{A_s} q_s(rs) d^2rs = 1,     ∬_{A_d} q_d(rd) d^2rd = 1.
%
%  Then the effective distance distribution is:
%     p_eff(rho) = ∬ ∬ q_s(rs) q_d(rd) δ(rho - |rs-rd|) d^2rs d^2rd
%
%  p_eff(rho) is a 1D PDF on [0, R_s + R_d] and satisfies ∫ p_eff(rho) d rho = 1.
%
%  ------------------------------------------------------------------------
%  Geometry/weight model used here
%  ------------------------------------------------------------------------
%  - Both source and detector are axisymmetric Gaussian weights:
%       w(r) = exp(-2 r^2 / w_mm^2)
%  - Each is truncated by a hard aperture:
%       0 <= r <= R_s   (source) ,  0 <= r <= R_d  (detector)
%  - We use the *normalized* 2D PDFs:
%       q(r) = w(r) / Z(R),   with Z(R) = ∬_{disk(R)} w(r) dA
%  - Sampling (method 1) uses inverse-CDF for the truncated 2D Gaussian.
%
%
%  Notes:
%    - This demo intentionally includes "naive" curves to show what *not* to do.
%    - For the repo, this file is meant as a reference implementation of p_eff(rho).
% ========================================================================

clc; clear; close all;

%% --------------------------- User settings ----------------------------
% Geometry / weights
w_mm = 0.50;     % Gaussian 1/e^2 radius in exp(-2 r^2 / w^2) [mm]
R_s  = 0.50;     % source aperture radius [mm]
R_d  = 0.50;     % detector aperture radius [mm]

% Monte-Carlo sampling for curve (1)
N_pairs_total = 3e6;   % sampling pairs for p_eff(rho)
N_plot        = 1e4;   % number of points to plot for geometry
Nbins         = 100;   % histogram bins for p_eff(rho)

% Baseline integral settings (method 2): radial quadrature resolution
NrS = 1100;
NrD = 1100;

% Hankel settings (method 3)
Nr_int = 4000;   % radial quadrature for Q(k)
Nk     = 2500;   % k grid points
kmax   = 350;    % [1/mm] upper k (increase if you want sharper agreement near rho=0)

% Numeric housekeeping
rng(1);          % reproducibility
if any([w_mm, R_s, R_d] <= 0), error('w_mm, R_s, R_d must be > 0'); end

%% --------------------- rho grid (shared by all curves) -----------------
rho_max  = R_s + R_d;
rho_edges = linspace(0, rho_max, Nbins+1);
drho      = rho_edges(2) - rho_edges(1);
rho_mid   = (rho_edges(1:end-1) + rho_edges(2:end))/2;
rho_mid   = rho_mid(:); % column

%% =======================================================================
%  0) Define the normalized 2D PDFs q_s(r), q_d(r) on the disks
% =======================================================================
% Unnormalized radial weight
wfun = @(r) exp(-2*(r.^2)/(w_mm^2));

% Normalization constant over a disk of radius R:
%   Z(R) = ∬ exp(-2 r^2 / w^2) dA = 2π ∫_0^R exp(-2 r^2 / w^2) r dr
%        = (π w^2 / 2) [1 - exp(-2 R^2 / w^2)]
Z = @(R) (pi*w_mm^2/2) * (1 - exp(-2*R.^2/w_mm^2));

Zs = Z(R_s);
Zd = Z(R_d);

% Normalized 2D densities (units 1/area) evaluated at radius r
qs = @(r) wfun(r) / Zs;   % valid for 0<=r<=R_s
qd = @(r) wfun(r) / Zd;   % valid for 0<=r<=R_d

% Radial PDFs for sampling radii (include Jacobian 2πr):
%   p_r(r) = 2π r q(r)
pr_s = @(r) 2*pi*r.*qs(r);
pr_d = @(r) 2*pi*r.*qd(r);

% Inverse CDF for the truncated 2D Gaussian radius:
%   CDF(r) = (1 - exp(-2 r^2/w^2)) / (1 - exp(-2 R^2/w^2))
inv_r_gauss_trunc = @(u, Rtr) (w_mm/sqrt(2)) * sqrt( -log( 1 - u*(1 - exp(-2*Rtr^2/w_mm^2)) ) );

%% =======================================================================
%  1) Sampled p_eff(rho)  (importance sampling from q_s and q_d)
% =======================================================================
u1  = rand(N_pairs_total,1);  th1 = 2*pi*rand(N_pairs_total,1);
u2  = rand(N_pairs_total,1);  th2 = 2*pi*rand(N_pairs_total,1);

rs  = inv_r_gauss_trunc(u1, R_s);
rd  = inv_r_gauss_trunc(u2, R_d);

xs  = rs .* cos(th1);  ys = rs .* sin(th1);
xd  = rd .* cos(th2);  yd = rd .* sin(th2);

rho_samp = hypot(xs - xd, ys - yd);

cnt    = histcounts(rho_samp, rho_edges);
p_samp = cnt ./ (sum(cnt)*drho + eps);
p_samp = p_samp(:);

% Store a small subset for geometry plot (left subplot)
src_plot = [xs(1:N_plot), ys(1:N_plot)];
det_plot = [xd(1:N_plot), yd(1:N_plot)];

%% =======================================================================
%  2) Baseline p_eff(rho) via real-space quadrature (manuscript definition)
% =======================================================================
% Strategy:
%   Sample radii r_s and r_d on grids, integrate over their radial PDFs.
%   For each (r_s, r_d), the relative angle φ is uniform on [0, 2π].
%   Distance is:
%       rho(φ) = sqrt(r_s^2 + r_d^2 - 2 r_s r_d cos φ)
%   For a given rho-bin [rho1, rho2], we can compute the *measure* of angles
%   φ that yield rho ∈ [rho1, rho2], without ever using a delta-function or
%   dealing with 1/sqrt singularities.
%
%   Then:
%     p_bin(b) = ∬ P( rho ∈ bin_b | r_s, r_d ) p_rS(r_s) p_rD(r_d) dr_s dr_d
%   and p_base(mid) = p_bin / drho.

rS = linspace(0, R_s, NrS).';   drS = rS(2)-rS(1);
rD = linspace(0, R_d, NrD).';   drD = rD(2)-rD(1);

% Radial PDFs on the grids (normalize numerically for robustness)
pRS = pr_s(rS);  pRS = pRS ./ max(trapz(rS, pRS), eps);
pRD = pr_d(rD);  pRD = pRD ./ max(trapz(rD, pRD), eps);

% Joint weights for the 2D radial quadrature (includes drS*drD)
Wrt = (pRS * pRD.') * (drS * drD);   % NrS x NrD

% Precompute some reusable matrices
rS_safe = rS; rS_safe(1) = eps;      % avoid division at r=0
rD_safe = rD; rD_safe(1) = eps;
[RR, TT] = ndgrid(rS_safe, rD_safe);
RT = 2*RR.*TT;

% Probability MASS in each rho bin
p_bin = zeros(numel(rho_edges)-1, 1);

for b = 1:numel(p_bin)
    rho1 = rho_edges(b);
    rho2 = rho_edges(b+1);

    % Feasible rho range for each (r,t):
    rmin = abs(RR - TT);
    rmax = RR + TT;

    % Clip requested interval to feasible interval
    a = max(rho1, rmin);
    c = min(rho2, rmax);

    ok = (c > a) & (RT > 0);   % non-empty interval and avoid r*t=0
    if ~any(ok(:))
        p_bin(b) = 0;
        continue;
    end

    % Map rho to cos(phi) via:
    %   rho^2 = r^2 + t^2 - 2 r t cos(phi)
    % => cos(phi) = (r^2 + t^2 - rho^2) / (2 r t)
    %
    % Note monotonicity: as rho increases, cos(phi) decreases.
    cos_a = (RR.^2 + TT.^2 - a.^2) ./ RT;   % corresponds to rho=a
    cos_c = (RR.^2 + TT.^2 - c.^2) ./ RT;   % corresponds to rho=c

    % Clamp to [-1,1] for safety
    cos_a = min(max(cos_a, -1), 1);
    cos_c = min(max(cos_c, -1), 1);

    % Angle measure on [0,2π] yielding rho in [a,c]:
    %   cos(phi) ∈ [cos_c, cos_a]
    % On [0,π], phi = acos(cos). There are two symmetric angles ±phi in [0,2π],
    % giving factor 2:
    dphi = 2 * (acos(cos_c) - acos(cos_a));
    dphi(~ok) = 0;

    % Conditional probability for this (r,t):
    prob = dphi / (2*pi);

    % Integrate over radii
    p_bin(b) = sum(sum(prob .* Wrt));
end

% Convert bin masses -> piecewise-constant PDF per bin
p_base = p_bin ./ (diff(rho_edges(:)) + eps);
p_base(~isfinite(p_base)) = 0;
p_base(p_base < 0) = 0;
p_base = p_base ./ (sum(p_base)*drho + eps);

%% =======================================================================
%  3) Hankel / Fourier–Bessel method (cross-check)
% =======================================================================
% For radial 2D PDFs q(r), the 2D Fourier transforms are:
%   Q(k) = ∬ q(r) e^{-i k·r} d^2r = 2π ∫_0^R q(r) J0(k r) r dr
% Then the displacement density is the convolution in 2D:
%   f_U(u) = (q_s * q_d)(u)  (in 2D),
% so in Fourier domain:  F_U(k) = Q_s(k) Q_d(k).
% The radial p_eff(rho) is then:
%   p_eff(rho) = ∫_{|u|=rho} f_U(u) dS
%              = rho ∫_0^∞ Q_s(k) Q_d(k) J0(k rho) k dk
%
% Numerically we truncate k to [0,kmax] and integrate with trapz.

k = linspace(0, kmax, Nk);

rS_h = linspace(0, R_s, Nr_int);
rD_h = linspace(0, R_d, Nr_int);

% Precompute q(r)*r for speed
qs_r = qs(rS_h) .* rS_h;
qd_r = qd(rD_h) .* rD_h;

Q_s = zeros(size(k));
Q_d = zeros(size(k));

for ii = 1:numel(k)
    ki = k(ii);
    Q_s(ii) = 2*pi * trapz(rS_h, qs_r .* besselj(0, ki*rS_h));
    Q_d(ii) = 2*pi * trapz(rD_h, qd_r .* besselj(0, ki*rD_h));
end

p_hankel = zeros(numel(rho_mid),1);
for j = 1:numel(rho_mid)
    r = rho_mid(j);
    integrand = (Q_s .* Q_d) .* besselj(0, k*r) .* k;
    p_hankel(j) = r * trapz(k, integrand);
end

p_hankel(~isfinite(p_hankel)) = 0;
p_hankel(p_hankel < 0) = 0;
p_hankel = p_hankel ./ (sum(p_hankel)*drho + eps);

%% =======================================================================
%  4) Two naive overlays (intentionally simplistic)
% =======================================================================
% These are common "quick guesses" people might try:
%   - start from uniform overlap disk distance PDF p_unif(rho) (equal radii)
%   - multiply by an ad hoc Gaussian decay in rho
%   - renormalize
%
% They are NOT generally correct, but are useful in the repo to illustrate
% why proper p_eff(rho) needs an actual derivation/integral.

Rref = min(R_s, R_d);
df   = 2*Rref;

p_unif = p_rho_circle(rho_mid, df);
p_unif = p_unif ./ (sum(p_unif)*drho + eps);

g1 = exp(-2*(rho_mid.^2)/(w_mm^2));
g2 = exp(-4*(rho_mid.^2)/(w_mm^2));

p_naive1 = p_unif .* g1;  p_naive1 = p_naive1 ./ (sum(p_naive1)*drho + eps);
p_naive2 = p_unif .* g2;  p_naive2 = p_naive2 ./ (sum(p_naive2)*drho + eps);

%% =======================================================================
%  5) Plot (repo-style: 2 subplots)
% =======================================================================
figure('Color','w');

% --- Left subplot: geometry points ---
subplot(1,2,1); hold on; box on; grid on; axis equal;

s1 = scatter(src_plot(:,1), src_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
s2 = scatter(det_plot(:,1), det_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

% Draw apertures for context
t = linspace(0,2*pi,512);
plot(R_s*cos(t), R_s*sin(t), 'k--', 'LineWidth', 1.0); % source aperture
plot(R_d*cos(t), R_d*sin(t), 'k-',  'LineWidth', 1.0); % detector aperture

xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Gaussian-weighted disks (overlap): w=%.3g mm, R_s=%.3g mm, R_d=%.3g mm', w_mm, R_s, R_d), ...
      'Interpreter','none');
legend([s1 s2], {'Input (source) points','Output (detector) points'}, 'Location','best');

lim = 1.10 * max(R_s, R_d);
xlim([-lim lim]);
ylim([-lim lim]);

% --- Right subplot: p_eff(rho) curves ---
subplot(1,2,2); hold on; box on; grid on;

% Histogram curve (sampled) as bars
bar(rho_mid, p_samp, 1.0, 'FaceAlpha', 0.35, 'EdgeColor','none');

% Overlays
plot(rho_mid, p_base,   'k-',  'LineWidth', 2.0);
plot(rho_mid, p_hankel, 'r--', 'LineWidth', 2.0);
plot(rho_mid, p_naive1, '-.',  'LineWidth', 1.8);
plot(rho_mid, p_naive2, ':',   'LineWidth', 1.8);

xlabel('\rho [mm]');
ylabel('p_{eff}(\rho)');
title('p_{eff}(\rho): sampled vs baseline vs Hankel (+ two naive overlays)');
xlim([0 rho_max]);

legend({ ...
    'Sampled (weighted pairs, hist)', ...
    'Baseline (real-space quadrature)', ...
    'Hankel (Fourier–Bessel cross-check)', ...
    'Naive: p_{unif}(\rho)\timesexp(-2\rho^2/w^2)', ...
    'Naive: p_{unif}(\rho)\timesexp(-4\rho^2/w^2)' ...
    }, 'Location','best');


%% =======================================================================
%  Local helper: uniform overlap disk distance PDF (Garcia-Pelayo)
%  df = diameter = 2R, support 0<=rho<=df
% =======================================================================
function p = p_rho_circle(rho, df)
    rho = double(rho(:));
    df  = double(df);

    p = zeros(size(rho));
    in = (rho >= 0) & (rho <= df);
    if ~any(in), return; end

    x = rho(in)./df;
    x = min(max(x, 0), 1);

    p(in) = (16.*rho(in))./(pi*df.^2).*acos(x) ...
          - (16.*rho(in).^2)./(pi*df.^3).*sqrt(max(0, 1 - x.^2));

    p(~isfinite(p)) = 0;
    p(p < 0) = 0;
end