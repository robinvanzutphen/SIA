%% ========================================================================
%  DEMO 3 — Annulus -> Annulus (concentric) uniform rings
%           p(rho) for distance between:
%             - one point uniformly sampled in an annulus (r_in..r_out)
%             - one point uniformly sampled in the SAME annulus (r_in..r_out)
%           with the SAME center (concentric)
%
%
%  Analytic method:
%    Covariogram / overlap-area identity:
%      p(rho) = (2*pi*rho)/(A_S*A_D) * Area( S ∩ (D shifted by rho) )
%    Here S and D are annuli. Use inclusion–exclusion with disk overlaps:
%      annulus = disk(r_out) \ disk(r_in)
%      A_cap(rho) = Ov(So,Do) - Ov(So,Di) - Ov(Si,Do) + Ov(Si,Di)
%
% ========================================================================

clc; clear; close all;

%% --------------------------- User settings ----------------------------
N_pairs_total = 1e7;   % total source-detector pairs to sample for p(rho)
N_plot        = 1e4;   % number of points to plot for geometry
Nbins         = 100;   % histogram bins for p(rho)

r_in  = 0.25;          % inner radius [mm]
r_out = 0.50;          % outer radius [mm]

rng(1); % reproducibility for the demo

%% --------------------- Histogram grid (rho) ----------------------------
rho_max  = 2*r_out;
edges    = linspace(0, rho_max, Nbins+1);
binw     = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

%% --------------------------- Main sampling ----------------------------
% Uniform sampling in an annulus:
%   theta ~ U(0,2pi)
%   r     ~ sqrt( (r_out^2 - r_in^2)*U + r_in^2 )

% --- Source points in annulus ---
th_s = 2*pi*rand(N_pairs_total,1);
r_s  = sqrt((r_out^2 - r_in^2)*rand(N_pairs_total,1) + r_in^2);
xs   = r_s .* cos(th_s);
ys   = r_s .* sin(th_s);

% --- Detector points in annulus ---
th_d = 2*pi*rand(N_pairs_total,1);
r_d  = sqrt((r_out^2 - r_in^2)*rand(N_pairs_total,1) + r_in^2);
xd   = r_d .* cos(th_d);
yd   = r_d .* sin(th_d);

% Pairwise distances rho
rho = hypot(xs - xd, ys - yd);

% Histogram counts and sampled density p(rho)
counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

%% -------------------- Analytic p(rho) on same grid ---------------------
pdf_ana = prho_analytic_annulus_annulus_cov(rho_cent, r_in, r_out, r_in, r_out);

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

% Annulus outlines for context
t = linspace(0, 2*pi, 512);
plot(r_out*cos(t), r_out*sin(t), 'k-',  'LineWidth', 1.0);
plot(r_in *cos(t), r_in *sin(t),  'k--', 'LineWidth', 1.0);

xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Annulus \\rightarrow Annulus (concentric): r_{in}=%.3g, r_{out}=%.3g mm', r_in, r_out), ...
      'Interpreter','tex');
legend([s1 s2], {'Input (annulus) points','Output (annulus) points'}, 'Location','best');

lim = 1.10*r_out;
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
function p = prho_analytic_annulus_annulus_cov(rho, Rs_in, Rs_out, Rd_in, Rd_out)
%PRHO_ANALYTIC_ANNULUS_ANNULUS_COV  Analytic p(rho) for annulus -> annulus (concentric)
%
% x ~ Unif{ Rs_in <= |x| <= Rs_out }
% y ~ Unif{ Rd_in <= |y| <= Rd_out }
% rho = |x - y|
%
% Covariogram identity:
%   p(rho) = (2*pi*rho)/(A_S*A_D) * A_cap(rho)
%
% Using inclusion–exclusion (annulus = outer disk minus inner disk):
%   A_cap = Ov(So,Do) - Ov(So,Di) - Ov(Si,Do) + Ov(Si,Di)

rho = double(rho(:));
p   = zeros(size(rho));

A_S = pi * (Rs_out^2 - Rs_in^2);
A_D = pi * (Rd_out^2 - Rd_in^2);

if A_S <= 0 || A_D <= 0
    error('Annulus areas must be positive. Check inner/outer radii.');
end

Acap = circle_overlap_area(rho, Rs_out, Rd_out) ...
     - circle_overlap_area(rho, Rs_out, Rd_in ) ...
     - circle_overlap_area(rho, Rs_in , Rd_out) ...
     + circle_overlap_area(rho, Rs_in , Rd_in );

p = (2*pi .* rho ./ (A_S * A_D)) .* Acap;

% Numerical safety
p(p < 0 & p > -1e-14) = 0;

p = reshape(p, size(rho));
end


function A = circle_overlap_area(d, R1, R2)
%CIRCLE_OVERLAP_AREA  Area of intersection of two disks at separation d
% (same helper as in DEMO 4; duplicated here for self-contained demo)

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

A(contained) = pi * min(R1, R2)^2;

dp = d(partial);
dp(dp == 0) = eps;

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