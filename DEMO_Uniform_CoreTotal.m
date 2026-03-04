%% ========================================================================
%  DEMO 2 — Core -> Total (concentric) uniform disks
%           p(rho) for distance between:
%             - one point uniformly sampled in a core disk (radius r_core)
%             - one point uniformly sampled in the full disk (radius R_total)
%           with the SAME center (overlapping, concentric)
%
%  Output:
%    - ONE figure with TWO subplots:
%        (Left)  sampled input/output points (plot 1e4 points, semi-transparent)
%        (Right) sampled p(rho) histogram (100 bins) + analytical p(rho)
%
%  Notes:
%    - We sample 1e7 pairs directly (no chunking).
%    - Support: 0 <= rho <= R_total + r_core
%
%  Function naming:
%    - p_rho_core_total(rho_grid, r_core, R_total)
% ========================================================================

clc; clear; close all;

%% --------------------------- User settings ----------------------------
N_pairs_total = 1e7;   % total source-detector pairs to sample for p(rho)
N_plot        = 1e4;   % number of points to plot for geometry
Nbins         = 100;   % histogram bins for p(rho)

R_total = 0.50;        % outer (total) radius [mm]
r_core  = 0.25;        % core radius [mm]

rng(1); % reproducibility for the demo

%% --------------------- Histogram grid (rho) ----------------------------
rho_max  = R_total + r_core;
edges    = linspace(0, rho_max, Nbins+1);
binw     = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

%% --------------------------- Main sampling ----------------------------
% Sample source points uniformly in CORE disk (radius r_core)
th_s = 2*pi*rand(N_pairs_total,1);
r_s  = r_core*sqrt(rand(N_pairs_total,1));
xs   = r_s .* cos(th_s);
ys   = r_s .* sin(th_s);

% Sample detector points uniformly in TOTAL disk (radius R_total)
th_d = 2*pi*rand(N_pairs_total,1);
r_d  = R_total*sqrt(rand(N_pairs_total,1));
xd   = r_d .* cos(th_d);
yd   = r_d .* sin(th_d);

% Pairwise distances rho
rho = hypot(xs - xd, ys - yd);

% Histogram counts and sampled density p(rho)
counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

%% --------------------------- Geometry points ---------------------------
src_plot = [xs(1:N_plot), ys(1:N_plot)];
det_plot = [xd(1:N_plot), yd(1:N_plot)];

%% ---------------------- Analytical p(rho) overlay ----------------------
pdf_theo = p_rho_core_total(rho_cent, r_core, R_total);
pdf_theo = pdf_theo(:).';
pdf_theo = pdf_theo ./ (sum(pdf_theo)*binw + eps); % robust normalization

%% ------------------------------- Plot ---------------------------------
figure('Color','w');

% --- Left subplot: geometry points ---
subplot(1,2,1); hold on; box on; grid on; axis equal;

s1 = scatter(src_plot(:,1), src_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
s2 = scatter(det_plot(:,1), det_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

% Draw outlines: total disk + core disk
t = linspace(0, 2*pi, 512);
plot(R_total*cos(t), R_total*sin(t), 'k-', 'LineWidth', 1.0);
plot(r_core *cos(t), r_core *sin(t),  'k--','LineWidth', 1.0);

xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Core \\rightarrow Total (concentric): r_{core}=%.3g, R=%.3g mm', r_core, R_total), ...
      'Interpreter','tex');
legend([s1 s2], {'Input (core) points','Output (total) points'}, 'Location','best');

lim = 1.10*R_total;
xlim([-lim lim]);
ylim([-lim lim]);

% --- Right subplot: sampled histogram + theoretical curve ---
subplot(1,2,2); hold on; box on; grid on;

bar(rho_cent, pdf_samp, 1.0, 'FaceAlpha', 0.35, 'EdgeColor','none');
plot(rho_cent, pdf_theo, 'k-', 'LineWidth', 2.0);

xlabel('\rho [mm]');
ylabel('p(\rho)');
title('Distance distribution: sampled (100 bins) vs analytical');
legend({'Sampled (MC)','Analytical'}, 'Location','best');
xlim([0 rho_max]);

%% ========================================================================
%  Local function: p_rho_core_total(rho_grid, r_core, R_total)
% ========================================================================
function p = p_rho_core_total(rho, r_core, R_total)
%P_RHO_CORE_TOTAL  Analytical distance PDF p(rho) for concentric disks:
%   - point A uniform in disk radius r_core
%   - point B uniform in disk radius R_total (R_total >= r_core)
%
% Support:
%   0 <= rho <= R_total + r_core
%
% Piecewise closed form (common in geometric probability):
%   For 0 <= rho <= (R_total - r_core):
%       p(rho) = 2*rho / R_total^2
%   For (R_total - r_core) < rho <= (R_total + r_core):
%       p(rho) = (rho/pi) * [ (2*alpha - sin(2*alpha))/R_total^2
%                           + (2*beta  - sin(2*beta ))/r_core^2 ]
%   where
%       alpha = acos( (r_core^2 + rho^2 - R_total^2) / (2*r_core*rho) )
%       beta  = acos( (R_total^2 + rho^2 - r_core^2) / (2*R_total*rho) )
%
% Numerical notes:
%   - handle rho=0 safely
%   - clamp acos arguments to [-1,1]

    rho_in = rho;               % keep original shape
    rho    = double(rho(:));
    a      = double(r_core);
    R      = double(R_total);

    p = zeros(size(rho));

    % Support mask
    in = (rho >= 0) & (rho <= (R + a));
    if ~any(in)
        p = reshape(p, size(rho_in));
        return;
    end

    r = rho(in);

    % Region 1: fully contained shifts
    r1 = (r > 0) & (r <= (R - a));   % exclude r=0 from division; set p(0)=0 anyway
    p_in = zeros(size(r));

    p_in(r1) = (2 .* r(r1)) ./ (R.^2);

    % Region 2: partial overlap
    r2 = (r > (R - a)) & (r <= (R + a)) & (r > 0);

    if any(r2)
        rr = r(r2);

        c1 = (a^2 + rr.^2 - R^2) ./ (2*a.*rr);
        c2 = (R^2 + rr.^2 - a^2) ./ (2*R.*rr);

        c1 = min(max(c1, -1), 1);
        c2 = min(max(c2, -1), 1);

        alpha = acos(c1);
        beta  = acos(c2);

        p_in(r2) = (rr/pi) .* ( (2*alpha - sin(2*alpha))./(R.^2) + (2*beta - sin(2*beta))./(a.^2) );
    end

    % Write back
    p(in) = p_in;

    % Safety cleanup
    p(~isfinite(p)) = 0;
    p(p < 0) = 0;

    % Return same shape as input
    p = reshape(p, size(rho_in));
end