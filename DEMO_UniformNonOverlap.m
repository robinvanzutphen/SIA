%% ========================================================================
%  DEMO 5 — Non-overlapping uniform disks (DRS-style)
%           p(rho) for the distance between:
%             - one point uniformly sampled in a source disk (radius rA)
%             - one point uniformly sampled in a detector disk (radius rB)
%           where the two disks are laterally separated by "sep" (no overlap).
%
%
%  Key idea for the semi-analytic p(rho):
%    The distance distribution can be written using a *covariogram / overlap* trick.
%
%    Let X ~ Uniform(disk A) and Y ~ Uniform(disk B). Define the displacement
%      U = X - Y.
%    Then the distance we care about is rho = |U|.
%
%    The PDF of U is the normalized overlap area between disk A and a shifted
%    disk B:
%      f_U(u) = Area( A ∩ (B + u) ) / (Area(A) * Area(B)).
%
%    This is because the convolution of uniform indicator functions of A and B
%    produces the overlap area.
%
%    To get the radial (distance) PDF p(rho), we integrate f_U(u) over the circle
%    of radius rho in displacement space:
%      p(rho) = ∫_{|u|=rho} f_U(u) dS
%            = rho * ∫_{0}^{2π} f_U( rho*[cosθ, sinθ] ) dθ.
%
%    When the disk centers are separated by "sep" along +x (B is shifted by [sep,0]),
%    the separation between the two disk centers inside the overlap term becomes:
%      s(θ) = | [sep,0] + rho*[cosθ, sinθ] |
%           = sqrt( sep^2 + rho^2 + 2*sep*rho*cosθ ).
%
%    Therefore:
%      p(rho) = rho / (Area(A)*Area(B)) * ∫_{0}^{2π} OverlapArea( s(θ); rA, rB ) dθ
%
%    Where OverlapArea(s; rA, rB) is the standard two-circle intersection area.
%
%  Notes:
%    - This is NOT a "closed-form" elementary expression for p(rho), but it is a
%      simple and fast 1D integral over θ (geometric, independent of MC transport).
%    - Support:
%        rho ∈ [ sep - (rA+rB),  sep + (rA+rB) ]
%      (for sep > rA+rB, the lower bound is > 0).
% ========================================================================

clc; clear; close all;

%% --------------------------- User settings ----------------------------
N_pairs_total = 1e7;   % total source-detector pairs to sample for p(rho)
N_plot        = 1e4;   % number of points to plot for geometry
Nbins         = 100;   % histogram bins for p(rho)

rA  = 0.30;            % source disk radius [mm]
rB  = 0.45;            % detector disk radius [mm]
sep = rA + rB + 0.25;  % center-to-center separation [mm] (choose > rA+rB for non-overlap)

% Numerical settings for the semi-analytic angular integral
Ntheta = 4096;         % angular samples for θ-integral (higher = smoother curve)

rng(1); % reproducibility for the demo

% Sanity check: enforce non-overlap
if sep <= (rA + rB)
    error('This demo is for NON-overlapping disks. Choose sep > rA + rB.');
end

%% --------------------- Histogram grid (rho) ----------------------------
rho_min = sep - (rA + rB);
rho_max = sep + (rA + rB);

edges    = linspace(rho_min, rho_max, Nbins+1);
binw     = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

%% --------------------------- Main sampling ----------------------------
% Disk A centered at (0,0); disk B centered at (sep,0)

% --- Sample points uniformly in disk A (source) ---
th_s = 2*pi*rand(N_pairs_total,1);
r_s  = rA*sqrt(rand(N_pairs_total,1));
xs   = r_s .* cos(th_s);
ys   = r_s .* sin(th_s);

% --- Sample points uniformly in disk B (detector) ---
th_d = 2*pi*rand(N_pairs_total,1);
r_d  = rB*sqrt(rand(N_pairs_total,1));
xd   = sep + r_d .* cos(th_d);    % shift by +sep in x
yd   =       r_d .* sin(th_d);

% Pairwise distances rho
rho = hypot(xs - xd, ys - yd);

% Histogram counts and sampled density p(rho)
counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

%% --------------------------- Geometry points ---------------------------
src_plot = [xs(1:N_plot), ys(1:N_plot)];
det_plot = [xd(1:N_plot), yd(1:N_plot)];

%% ---------------------- Semi-analytic p(rho) overlay -------------------
pdf_theo = p_rho_nonoverlap_disks(rho_cent, sep, rA, rB, Ntheta);

% Robust normalization (should already be very close)
pdf_theo = pdf_theo(:).';
pdf_theo = pdf_theo ./ (sum(pdf_theo)*binw + eps);

%% ------------------------------- Plot ---------------------------------
figure('Color','w');

% --- Left subplot: geometry points ---
subplot(1,2,1); hold on; box on; grid on; axis equal;

s1 = scatter(src_plot(:,1), src_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
s2 = scatter(det_plot(:,1), det_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

% Draw disk outlines for context
t = linspace(0, 2*pi, 512);
plot(rA*cos(t),          rA*sin(t),          'k--', 'LineWidth', 1.0);  % disk A
plot(sep + rB*cos(t),    rB*sin(t),          'k-',  'LineWidth', 1.0);  % disk B

xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Non-overlapping disks: r_A=%.3g, r_B=%.3g, sep=%.3g mm', rA, rB, sep), ...
      'Interpreter','tex');
legend([s1 s2], {'Input (disk A) points','Output (disk B) points'}, 'Location','best');

limx = [min(-rA, sep-rB) , max(rA, sep+rB)];
pad  = 0.15*(limx(2)-limx(1));
xlim([limx(1)-pad, limx(2)+pad]);
ylim([-max(rA,rB)-pad, max(rA,rB)+pad]);

% --- Right subplot: sampled histogram + semi-analytic curve ---
subplot(1,2,2); hold on; box on; grid on;

bar(rho_cent, pdf_samp, 1.0, 'FaceAlpha', 0.35, 'EdgeColor','none');
plot(rho_cent, pdf_theo, 'k-', 'LineWidth', 2.0);

xlabel('\rho [mm]');
ylabel('p(\rho)');
title('Distance distribution: sampled (100 bins) vs semi-analytic (1D integral)');
legend({'Sampled (MC)','Semi-analytic (covariogram + overlap)'}, 'Location','best');
xlim([rho_min rho_max]);

%% ========================================================================
%  Local functions
% ========================================================================

function p = p_rho_nonoverlap_disks(rho, sep, rA, rB, Ntheta)
%P_RHO_NONOVERLAP_DISKS  Semi-analytic p(rho) for two uniform disks separated by sep.
%
%   p(rho) = rho/(AreaA*AreaB) * ∫_0^{2π} OverlapArea(s(θ); rA, rB) dθ
%   with s(θ) = sqrt(sep^2 + rho^2 + 2*sep*rho*cosθ).
%
% Inputs:
%   rho   : vector of distances
%   sep   : center-to-center separation (disk B is at [sep,0], disk A at [0,0])
%   rA    : radius of disk A
%   rB    : radius of disk B
%   Ntheta: number of θ samples for the 1D integral
%
% Output:
%   p     : p(rho) evaluated at rho (same size)

    rho_in = rho;
    rho    = double(rho(:));
    sep    = double(sep);
    rA     = double(rA);
    rB     = double(rB);

    AreaA = pi*rA^2;
    AreaB = pi*rB^2;

    % θ grid for the 1D integral
    theta = linspace(0, 2*pi, Ntheta);
    cth   = cos(theta);

    p = zeros(size(rho));

    % Support (works for all sep, but for non-overlap the lower bound is >0)
    rho_min = abs(sep - (rA + rB));
    rho_max = sep + (rA + rB);

    in = (rho >= rho_min) & (rho <= rho_max);
    if ~any(in)
        p = reshape(p, size(rho_in));
        return;
    end

    % Evaluate integral for each rho (loop over rho is fine for 100-bin demo)
    for i = find(in).'
        r = rho(i);

        % center separation inside overlap area as a function of θ
        s = sqrt(sep^2 + r^2 + 2*sep*r*cth);

        % overlap area between two disks with radii rA and rB separated by s
        Aov = overlap_area_two_circles(s, rA, rB);

        % angular integral
        I = trapz(theta, Aov);

        % radial PDF p(rho)
        p(i) = (r / (AreaA*AreaB)) * I;
    end

    % Safety cleanup
    p(~isfinite(p)) = 0;
    p(p < 0) = 0;

    % Return in original shape
    p = reshape(p, size(rho_in));
end

function A = overlap_area_two_circles(s, a, b)
%OVERLAP_AREA_TWO_CIRCLES  Area of intersection of two circles (radii a, b)
% separated by distance s (vectorized over s).
%
% Standard piecewise result:
%   - No overlap: s >= a+b   => A=0
%   - Full containment: s <= |a-b| => A=pi*min(a,b)^2
%   - Partial overlap: use acos terms - 0.5*sqrt(…) term

    s = double(s(:));
    a = double(a);
    b = double(b);

    A = zeros(size(s));

    % No overlap
    no = (s >= (a + b));
    A(no) = 0;

    % One contains the other
    co = (s <= abs(a - b));
    A(co) = pi * min(a,b)^2;

    % Partial overlap
    pa = ~(no | co);
    sp = s(pa);

    % acos arguments with clamping
    c1 = (sp.^2 + a^2 - b^2) ./ (2*sp*a);
    c2 = (sp.^2 + b^2 - a^2) ./ (2*sp*b);
    c1 = min(max(c1, -1), 1);
    c2 = min(max(c2, -1), 1);

    term = (-sp + a + b) .* (sp + a - b) .* (sp - a + b) .* (sp + a + b);
    term = max(term, 0); % numeric safety

    A(pa) = a^2 .* acos(c1) + b^2 .* acos(c2) - 0.5 .* sqrt(term);

    % reshape back to match input shape (vector column is fine for our use)
end