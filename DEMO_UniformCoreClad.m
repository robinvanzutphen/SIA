%% ========================================================================
%  DEMO 4 — Core -> Cladding (concentric) uniform disk -> uniform annulus
%           p(rho) for distance between:
%             - one point uniformly sampled in the core disk   (radius r_core)
%             - one point uniformly sampled in the cladding annulus (r_core..R_out)
%           with the SAME center (concentric)
%
%
%  Analytic method:
%    Covariogram / overlap-area identity:
%      p(rho) = (2*pi*rho)/(A_S*A_D) * Area( S ∩ (D shifted by rho) )
%    Here:
%      S = disk(R_core)
%      D = annulus(R_core..R_out) = disk(R_out) \ disk(R_core)
%    => overlap area:
%      A_cap(rho) = Ov(disk(R_core), disk(R_out)) - Ov(disk(R_core), disk(R_core))
%
% ========================================================================

clc; clear; close all;

%% --------------------------- User settings ----------------------------
N_pairs_total = 1e7;   % total source-detector pairs to sample for p(rho)
N_plot        = 1e4;   % number of points to plot for geometry
Nbins         = 100;   % histogram bins for p(rho)

r_core = 0.25;         % core radius [mm]
R_out  = 0.50;         % total (outer) radius [mm]  => cladding is [r_core, R_out]

rng(1); % reproducibility for the demo

%% --------------------- Histogram grid (rho) ----------------------------
rho_max  = r_core + R_out;
edges    = linspace(0, rho_max, Nbins+1);
binw     = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

%% --------------------------- Main sampling ----------------------------
% Uniform sampling in a disk:
%   theta ~ U(0,2pi)
%   r     ~ R*sqrt(U)
%
% Uniform sampling in an annulus (r_in..r_out):
%   theta ~ U(0,2pi)
%   r     ~ sqrt( (r_out^2 - r_in^2)*U + r_in^2 )

% --- Source points in core disk ---
th_s = 2*pi*rand(N_pairs_total,1);
r_s  = r_core * sqrt(rand(N_pairs_total,1));
xs   = r_s .* cos(th_s);
ys   = r_s .* sin(th_s);

% --- Detector points in cladding annulus [r_core, R_out] ---
th_d = 2*pi*rand(N_pairs_total,1);
r_d  = sqrt((R_out^2 - r_core^2)*rand(N_pairs_total,1) + r_core^2);
xd   = r_d .* cos(th_d);
yd   = r_d .* sin(th_d);

% Pairwise distances rho
rho = hypot(xs - xd, ys - yd);

% Histogram counts and sampled density p(rho)
counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

%% -------------------- Analytic p(rho) on same grid ---------------------
pdf_ana = prho_analytic_disk_annulus_cov(rho_cent, r_core, r_core, R_out);

%% --------------------------- Geometry points ---------------------------
src_plot = [xs(1:N_plot), ys(1:N_plot)];
det_plot = [xd(1:N_plot), yd(1:N_plot)];

%% ------------------------------- Plot ---------------------------------
figure('Color','w');

% --- Left subplot: geometry points ---
subplot(1,2,1); hold on; box on; grid on; axis equal;

s1 = scatter(src_plot(:,1), src_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
s2 = scatter(det_plot(:,1), det_plot(:,2), 6, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

% Draw outlines: outer radius + core radius (cladding is the ring between them)
t = linspace(0, 2*pi, 512);
plot(R_out*cos(t),  R_out*sin(t),  'k-',  'LineWidth', 1.0);
plot(r_core*cos(t), r_core*sin(t), 'k--', 'LineWidth', 1.0);

xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Core \\rightarrow Cladding (concentric): r_{core}=%.3g, R_{out}=%.3g mm', r_core, R_out), ...
      'Interpreter','tex');
legend([s1 s2], {'Input (core disk) points','Output (cladding annulus) points'}, 'Location','best');

lim = 1.10*R_out;
xlim([-lim lim]);
ylim([-lim lim]);

% --- Right subplot: sampled histogram + analytic overlay ---
subplot(1,2,2); hold on; box on; grid on;

bar(rho_cent, pdf_samp, 1.0, 'FaceAlpha', 0.35, 'EdgeColor','none');
plot(rho_cent, pdf_ana, 'k-', 'LineWidth', 2.0);

xlabel('\rho [mm]');
ylabel('p(\rho)');
title('Distance distribution: sampled (bars) + analytic (line)');
legend({'Sampled (MC)','Analytic (covariogram)'}, 'Location','best');
xlim([0 rho_max]);

%% ============================ Functions ================================
function p = prho_analytic_disk_annulus_cov(rho, R_disk, R_ann_in, R_ann_out)
%PRHO_ANALYTIC_DISK_ANNULUS_COV  Analytic p(rho) for disk -> annulus (concentric)
%
%   x ~ Unif{ |x| <= R_disk }
%   y ~ Unif{ R_ann_in <= |y| <= R_ann_out }
%   rho = |x - y|
%
% Covariogram identity:
%   p(rho) = (2*pi*rho)/(A_S*A_D) * A_cap(rho),
% where A_cap(rho) is the overlap area between S and D shifted by distance rho.
%
% Annulus = disk(R_ann_out) \ disk(R_ann_in):
%   A_cap(rho) = Ov(disk(R_disk), disk(R_ann_out)) - Ov(disk(R_disk), disk(R_ann_in))

rho = double(rho(:));
p   = zeros(size(rho));

A_S = pi * R_disk^2;
A_D = pi * (R_ann_out^2 - R_ann_in^2);

if A_S <= 0 || A_D <= 0
    error('Areas must be positive. Check radii.');
end

Acap = circle_overlap_area(rho, R_disk, R_ann_out) ...
     - circle_overlap_area(rho, R_disk, R_ann_in);

p = (2*pi .* rho ./ (A_S * A_D)) .* Acap;

% Numerical safety: clip tiny negatives due to floating-point rounding
p(p < 0 & p > -1e-14) = 0;

p = reshape(p, size(rho));
end


function A = circle_overlap_area(d, R1, R2)
%CIRCLE_OVERLAP_AREA  Area of intersection of two disks at separation d
%
% Vectorized piecewise definition:
%   - d >= R1+R2   : 0
%   - d <= |R1-R2| : area of smaller disk
%   - else         : lens formula (acos + sqrt)
%
% Robustness:
%   - clamps acos arguments to [-1,1]
%   - guards against negative sqrt due to rounding

d  = double(d(:));
R1 = double(R1);
R2 = double(R2);

A = zeros(size(d));

if R1 < 0 || R2 < 0
    error('Radii must be non-negative.');
end
if R1 == 0 || R2 == 0
    return;
end

noOverlap = d >= (R1 + R2);
contained = d <= abs(R1 - R2);
partial   = ~(noOverlap | contained);

% Containment: smaller disk fully inside larger
A(contained) = pi * min(R1, R2)^2;

% Partial overlap (lens)
dp = d(partial);
dp(dp == 0) = eps; % safety (should not occur in partial regime)

arg1 = (dp.^2 + R1^2 - R2^2) ./ (2*dp*R1);
arg2 = (dp.^2 + R2^2 - R1^2) ./ (2*dp*R2);
arg1 = max(-1, min(1, arg1));
arg2 = max(-1, min(1, arg2));

term1 = R1^2 .* acos(arg1);
term2 = R2^2 .* acos(arg2);

rootTerm = (-dp + R1 + R2) .* (dp + R1 - R2) .* (dp - R1 + R2) .* (dp + R1 + R2);
rootTerm = max(0, rootTerm);

term3 = 0.5 .* sqrt(rootTerm);

A(partial) = term1 + term2 - term3;

A = reshape(A, size(d));
end