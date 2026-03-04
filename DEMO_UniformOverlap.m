%% ========================================================================
%  DEMO 1 — Overlapping (coaxial) uniform disk source + uniform disk detector
%           p(rho) for two random points drawn uniformly from the SAME disk
% ========================================================================

clc; clear; close all;

%% --------------------------- User settings ----------------------------
N_pairs_total = 1e7;   % total source-detector pairs to sample for p(rho)
N_plot        = 1e4;   % number of points to plot for geometry
Nbins         = 100;   % histogram bins for p(rho)

diam = 1.0;            % disk diameter (units arbitrary but consistent; e.g., mm)
R    = diam/2;         % disk radius

rng(1); % reproducibility for the demo

%% --------------------------- Distance histogram ------------------------
% Histogram edges and bin width: support is rho in [0, diam]
edges  = linspace(0, diam, Nbins+1);
binw   = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

% Sample two *independent* uniform points in the same disk.
th_s = 2*pi*rand(N_pairs_total,1);
r_s  = R*sqrt(rand(N_pairs_total,1));
xs   = r_s .* cos(th_s);
ys   = r_s .* sin(th_s);

th_d = 2*pi*rand(N_pairs_total,1);
r_d  = R*sqrt(rand(N_pairs_total,1));
xd   = r_d .* cos(th_d);
yd   = r_d .* sin(th_d);

% Pairwise distances rho
rho = hypot(xs - xd, ys - yd);

% Histogram counts and sampled density p(rho)
counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

%% --------------------------- Geometry points ---------------------------
% For the left subplot we only plot a subset (N_plot) to keep it readable.
% We reuse the already sampled coordinates (first N_plot samples).
src_plot = [xs(1:N_plot), ys(1:N_plot)];
det_plot = [xd(1:N_plot), yd(1:N_plot)];

%% ---------------------- Theoretical p(rho) overlay ---------------------
pdf_theo = p_rho_uniform_overlap(rho_cent, diam);
pdf_theo = pdf_theo(:).';
pdf_theo = pdf_theo ./ (sum(pdf_theo)*binw + eps); % robust normalization

%% ------------------------------- Plot ---------------------------------
figure('Color','w');

% --- Left subplot: geometry points ---
subplot(1,2,1); hold on; box on; grid on; axis equal;

% Small markers + semi-transparent look:
% (If MarkerFaceAlpha is not supported in your MATLAB version, remove it.)
s1 = scatter(src_plot(:,1), src_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
s2 = scatter(det_plot(:,1), det_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

% Disk outline for context
t = linspace(0, 2*pi, 512);
plot(R*cos(t), R*sin(t), 'k-', 'LineWidth', 1.0);

xlabel('x');
ylabel('y');
title(sprintf('Uniform overlap disk: diam = %.3g', diam), 'Interpreter','none');
legend([s1 s2], {'Input (source) points','Output (detector) points'}, 'Location','best');

lim = 1.1*R;
xlim([-lim lim]);
ylim([-lim lim]);

% --- Right subplot: sampled histogram + theoretical curve ---
subplot(1,2,2); hold on; box on; grid on;

bar(rho_cent, pdf_samp, 1.0, 'FaceAlpha', 0.35, 'EdgeColor','none');
plot(rho_cent, pdf_theo, 'k-', 'LineWidth', 2.0);

xlabel('\rho');
ylabel('p(\rho)');
title('Distance distribution: sampled (100 bins) vs theoretical');
legend({'Sampled (MC)','Theoretical'}, 'Location','best');
xlim([0 diam]);

%% ========================================================================
%  Local function: p_rho_uniform_overlap(rho_grid, diam)
% ========================================================================
function p = p_rho_uniform_overlap(rho, diam)
%P_RHO_UNIFORM_OVERLAP  Distance PDF p(rho) for two independent uniform points
%                       in the SAME disk of diameter "diam".
%
%   p = p_rho_uniform_overlap(rho_grid, diam)
%
% Inputs
%   rho  : vector of distances at which to evaluate p(rho)
%   diam : disk diameter (same units as rho)
%
% Output
%   p    : p(rho) evaluated at rho (same size as rho)
%
% Closed-form result (support 0 <= rho <= diam):
%   p(rho) = (16*rho)/(pi*diam^2) * acos(rho/diam)
%            - (16*rho^2)/(pi*diam^3) * sqrt(1 - (rho/diam)^2)

    rho_in = rho;                 % keep original shape
    rho    = double(rho(:));
    d      = double(diam);

    p = zeros(size(rho));

    in = (rho >= 0) & (rho <= d);
    if any(in)
        x = rho(in) ./ d;
        x = min(max(x, 0), 1); % numeric safety

        p(in) = (16.*rho(in))./(pi*d.^2).*acos(x) ...
              - (16.*rho(in).^2)./(pi*d.^3).*sqrt(max(0, 1 - x.^2));
    end

    % Safety cleanup
    p(~isfinite(p)) = 0;
    p(p < 0) = 0;

    % Return same shape as input rho
    p = reshape(p, size(rho_in));
end