%% ========================================================================
%  DEMO 6 — Square -> Square (touching corners, top-right)
%           p(rho) for distance between:
%             - one point uniformly sampled in a source square (side s_src)
%             - one point uniformly sampled in a detector square (side s_det)
%           where the detector square is shifted by [dx, dy] relative to the
%           source (centers offset).
%
%  This demo uses the same "covariogram / overlap" trick as the non-overlapping
%  disk case, but now the overlap area is *closed form* and extremely cheap.
%
%
%  Let X ~ Uniform(source square S) and Y ~ Uniform(detector square D).
%  Define displacement U = X - Y. Then rho = |U|.
%
%  The 2D PDF of U is proportional to the overlap area between S and D shifted
%  by u:
%      f_U(u) = Area( S ∩ (D + u) ) / (Area(S) * Area(D)).
%
%  To get the radial distance PDF p(rho), integrate f_U over the circle of
%  radius rho in displacement space:
%      p(rho) = ∫_{|u|=rho} f_U(u) dS
%            = rho * ∫_0^{2π} f_U( rho*[cosθ, sinθ] ) dθ.
%
%  If D is centered at offset d = [dx, dy] relative to S, then inside the
%  overlap we evaluate:
%      u(θ) = d + rho*[cosθ, sinθ].
%
%  For axis-aligned squares, Area( S ∩ (S + t) ) is separable:
%      A_cap(t) = max(0, s - |t_x|) * max(0, s - |t_y|)
%  (and similarly for unequal square sizes, using the 1D interval overlap).
%
%  Therefore p(rho) becomes a simple 1D integral over θ of a cheap expression.
%
%  Notes:
%    - This is a *semi-analytic* 1D θ-integral (not "elementary closed-form"),
%      but it is geometrically exact and computationally negligible vs MC.
% ========================================================================

clc; clear; close all;

%% --------------------------- User settings ----------------------------
N_pairs_total = 1e7;   % total source-detector pairs to sample for p(rho)
N_plot        = 1e4;   % number of points to plot for geometry
Nbins         = 100;   % histogram bins for p(rho)

% --- Geometry (your case: equal squares touching corners) ---
s_src = 0.50;          % source square side [mm]
s_det = 0.50;          % detector square side [mm]
dx    = s_src;         % detector center offset in x [mm] (top-right touching corners)
dy    = s_src;         % detector center offset in y [mm]

% --- Numerical settings for semi-analytic θ integral ---
Ntheta = 4096;         % angular integration resolution (higher = smoother)

rng(1); % reproducibility for the demo

%% --------------------- Histogram grid (rho) ----------------------------
% Safe support upper bound: farthest corner-to-corner distance
rho_max  = hypot(dx + (s_src/2 + s_det/2), dy + (s_src/2 + s_det/2));
edges    = linspace(0, rho_max, Nbins+1);
binw     = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

%% --------------------------- Main sampling ----------------------------
% Source square centered at (0,0)
Ps = sample_square(N_pairs_total, s_src);

% Detector square centered at (dx,dy)
Pd = sample_square(N_pairs_total, s_det) + [dx, dy];

% Pairwise distances rho
rho = hypot(Ps(:,1) - Pd(:,1), Ps(:,2) - Pd(:,2));

% Histogram counts and sampled density p(rho)
counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

%% --------------------------- Geometry points ---------------------------
src_plot = Ps(1:N_plot, :);
det_plot = Pd(1:N_plot, :);

%% ---------------------- Semi-analytic p(rho) overlay -------------------
pdf_theo = p_rho_square_square_cov(rho_cent, s_src, s_det, dx, dy, Ntheta);

% Robust normalization (finite θ discretization can introduce small drift)
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

% Draw square outlines for context
draw_square_outline(0, 0, s_src, 'k--');       % source
draw_square_outline(dx, dy, s_det, 'k-');      % detector

xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Square \\rightarrow Square: s_{src}=%.3g, s_{det}=%.3g mm, offset=[%.3g %.3g] mm', ...
      s_src, s_det, dx, dy), 'Interpreter','tex');
legend([s1 s2], {'Input (source) points','Output (detector) points'}, 'Location','best');

% Reasonable axis limits
xmin = min(-s_src/2, dx - s_det/2);
xmax = max( s_src/2, dx + s_det/2);
ymin = min(-s_src/2, dy - s_det/2);
ymax = max( s_src/2, dy + s_det/2);
pad  = 0.15 * max([xmax-xmin, ymax-ymin]);
xlim([xmin-pad, xmax+pad]);
ylim([ymin-pad, ymax+pad]);

% --- Right subplot: sampled histogram + semi-analytic curve ---
subplot(1,2,2); hold on; box on; grid on;

bar(rho_cent, pdf_samp, 1.0, 'FaceAlpha', 0.35, 'EdgeColor','none');
plot(rho_cent, pdf_theo, 'k-', 'LineWidth', 2.0);

xlabel('\rho [mm]');
ylabel('p(\rho)');
title('Distance distribution: sampled (100 bins) vs semi-analytic (covariogram, 1D \theta integral)');
legend({'Sampled (MC)','Semi-analytic (covariogram)'}, 'Location','best');
xlim([0 rho_max]);

%% ========================================================================
%  Local functions
% ========================================================================

function P = sample_square(N, side_mm)
%SAMPLE_SQUARE  Uniform points in a square of side "side_mm" centered at origin.
    half = side_mm/2;
    P = (rand(N,2) - 0.5) * (2*half);
end

function draw_square_outline(cx, cy, side_mm, ls)
%DRAW_SQUARE_OUTLINE  Draw outline of an axis-aligned square centered at (cx,cy).
    h = side_mm/2;
    x = cx + [-h, +h, +h, -h, -h];
    y = cy + [-h, -h, +h, +h, -h];
    plot(x, y, ls, 'LineWidth', 1.0);
end

function p = p_rho_square_square_cov(rho, s_src, s_det, dx, dy, Ntheta)
%P_RHO_SQUARE_SQUARE_COV  Semi-analytic p(rho) via covariogram + 1D θ integral.
%
% General (possibly unequal sizes) axis-aligned squares:
%   Source square S: side s_src, centered at (0,0)
%   Detector square D: side s_det, centered at (dx,dy)
%
% p(rho) = rho/(A_S*A_D) * ∫_0^{2π} A_cap( t_x(θ), t_y(θ) ) dθ
% where
%   [t_x(θ), t_y(θ)] = [dx,dy] + rho*[cosθ, sinθ]
% and A_cap is overlap area between two squares with relative shift (t_x,t_y).
%
% For axis-aligned squares, overlap area factorizes into independent 1D overlaps:
%   overlap_x = max(0, (s_src/2 + s_det/2) - |t_x|)
%   overlap_y = max(0, (s_src/2 + s_det/2) - |t_y|)
%   A_cap = overlap_x * overlap_y

    rho_in = rho;
    rho    = double(rho(:));

    hs = double(s_src)/2;
    hd = double(s_det)/2;

    % Areas
    As = double(s_src)^2;
    Ad = double(s_det)^2;

    % θ grid for the 1D integral
    theta = linspace(0, 2*pi, Ntheta);
    ct = cos(theta);
    st = sin(theta);

    % Build [Nrho x Ntheta] arrays via implicit expansion (R2016b+)
    tx = double(dx) + rho * ct;
    ty = double(dy) + rho * st;

    % 1D interval overlap lengths
    ovx = max(0, (hs + hd) - abs(tx));
    ovy = max(0, (hs + hd) - abs(ty));

    % 2D overlap area (covariogram)
    Acap = ovx .* ovy;

    % θ integral for each rho
    Itheta = trapz(theta, Acap, 2);

    % Radial PDF
    p = (rho ./ (As*Ad)) .* Itheta;

    % Safety cleanup
    p(~isfinite(p)) = 0;
    p(p < 0) = 0;

    % Return same shape
    p = reshape(p, size(rho_in));
end