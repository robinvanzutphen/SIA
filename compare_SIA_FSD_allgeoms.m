

clc; clear; close all;

%% ========================================================================
%  0) CASES TO ANALYZE (edit here)
% ========================================================================
% Index map (default order below):
%  1) disk -> disk (overlap, equal)
%  2) core -> total (concentric)
%  3) annulus -> annulus (concentric)
%  4) core -> cladding (concentric)
%  5) non-overlapping disks
%  6) square -> square (touching corners, top-right)
%  7) Gaussian source + Gaussian-weighted detector (overlap, cap-fix)
CasesToRun = [1 1 1 1 1 1 1];   % 1=run, 0=skip

%% ========================================================================
%  1) Global settings (shared by ALL cases unless explicitly noted)
% ========================================================================

% ---- Optics ----
NA    = 0.22;
n_med = 1.35;
theta_acc_max = asin(NA/n_med);

mu_a  = 0.10;              % [1/mm]
mu_sp = 5.00;              % [1/mm]
g     = 0.90;
mu_s  = mu_sp/(1-g);

% ---- Monte Carlo ----
reps_per_cfg = 2;          % set 1 for quick; >=5 for nicer distributions
nphoton      = 2e6;
gpu_id       = '11';

% ---- Domain ----
mm_per_voxel = 1;          % 1 voxel = 1 mm
Lxy = 240; Lz = 160;
vol = uint8(ones(Lxy,Lxy,Lz));

src_xy = round(Lxy/2);
src_z  = 0;
srcpos = [src_xy src_xy src_z];
srcdir = [0 0 1];
src_center_mm = double(srcpos(1:2)) * mm_per_voxel;

% ---- Pathlength bins (shape-only PDFs) ----
L_max  = 5;      % [mm]
nLbins = 1000;
edges_L = linspace(0, L_max, nLbins+1);
L_mid   = (edges_L(1:end-1)+edges_L(2:end))/2;
dL      = edges_L(2)-edges_L(1);

% ---- rho bins for mixing and SIA p(L) ----
nbins_rho = 5000;

% ---- Geometry sampling for plots + sanity p(rho) ----
Nbins_prho_plot   = 50;      % Fig 2 sampled curve: fixed 50 bins
Ngeom_plot        = 5000;    % store this many points for Fig 1
Ngeom_pairs_prho  = 2e6;     % pairs for sampled p(rho) sanity curve (Fig 2)
Ngeom_pairs_mix   = 1e7;     % ONLY when p(rho)=sampled-only (annulus cases + square)

%% ========================================================================
%  2) MCX base config (shared)
% ========================================================================

cfg0 = struct();
cfg0.nphoton      = nphoton;
cfg0.vol          = vol;
cfg0.unitinmm     = mm_per_voxel;

cfg0.srcpos       = srcpos;
cfg0.srcdir       = srcdir;

cfg0.autopilot    = 1;
cfg0.tstart       = 0;
cfg0.tend         = 5e-9;
cfg0.tstep        = 5e-9;

cfg0.isreflect    = 0;
cfg0.respin       = 1;
cfg0.issave2pt    = 0;
cfg0.savedetflag  = 'dspvx';
cfg0.maxdetphoton = nphoton;

cfg0.gpuid        = gpu_id;
cfg0.issrcfrom0   = 1;
cfg0.seed         = -1;

% Finite NA launch (encoded in photon transport)
cfg0.angleinvcdf  = linspace(0, theta_acc_max/pi, 5);

% Boundary conditions (NO SPACES) — use -z face as detector plane
cfg0.bc           = 'aa_aaa001000';

% Medium
cfg0.prop         = [0 0 1 1; mu_a mu_s g n_med];

Nlaunch = double(cfg0.nphoton * cfg0.respin);

%% ========================================================================
%  3) Define cases (1..7)
% ========================================================================

% --- Disk radii for script-1 cases ---
R_det  = 0.50;   % outer radius [mm]
r_core = 0.25;   % core radius [mm]

% --- Non-overlapping disks ---
rA  = 0.30;                       % source radius [mm]
rB  = 0.45;                       % detector radius [mm]
sep = rA + rB + 0.25*R_det;       % [mm] > rA+rB

% --- Square case ---
side_src = 0.50;  % [mm]
side_det = 0.50;  % [mm]
det_offset_sq_mm = [side_src, side_src]; % top-right touching corners

% --- Gaussian case (non-uniform) ---
w_mm     = 0.50;   % 1/e^2 radius in exp(-2 r^2 / w^2) [mm]
R_det_g  = 0.50;   % detector disk aperture radius [mm]
Nw_cap   = 10;     % source cap multiplier: R_src_cap = R_det + Nw_cap*w

CasesAll = { ...
  struct('name','1) disk -> disk (overlap, equal)', ...
         'type','uniform_shapes', ...
         'src_shape','disk',    'src_rin',0,      'src_rout',R_det, 'src_side',[], ...
         'det_shape','disk',    'det_rin',0,      'det_rout',R_det, 'det_side',[], ...
         'det_offset_mm',[0 0], ...
         'p_mode','analytic_circle'), ...

  struct('name','2) core -> total (concentric)', ...
         'type','uniform_shapes', ...
         'src_shape','disk',    'src_rin',0,      'src_rout',r_core, ...
         'det_shape','disk',    'det_rin',0,      'det_rout',R_det, ...
         'det_offset_mm',[0 0], ...
         'p_mode','analytic_core_total'), ...

  struct('name','3) annulus -> annulus (concentric)', ...
         'type','uniform_shapes', ...
         'src_shape','annulus', 'src_rin',r_core, 'src_rout',R_det, ...
         'det_shape','annulus', 'det_rin',r_core, 'det_rout',R_det, ...
         'det_offset_mm',[0 0], ...
         'p_mode','sampled_only'), ...

  struct('name','4) core -> cladding (concentric)', ...
         'type','uniform_shapes', ...
         'src_shape','disk',    'src_rin',0,      'src_rout',r_core, ...
         'det_shape','annulus', 'det_rin',r_core, 'det_rout',R_det, ...
         'det_offset_mm',[0 0], ...
         'p_mode','sampled_only'), ...

  struct('name','5) non-overlapping disks', ...
         'type','uniform_shapes', ...
         'src_shape','disk',    'src_rin',0,      'src_rout',rA, ...
         'det_shape','disk',    'det_rin',0,      'det_rout',rB, ...
         'det_offset_mm',[sep 0], ...
         'p_mode','analytic_DRS'), ...

  struct('name','6) square 0.5mm -> square 0.5mm (touching corners, top-right)', ...
         'type','square_square', ...
         'src_shape','square',  'src_rin',[], 'src_rout',[], 'src_side',side_src, ...
         'det_shape','square',  'det_rin',[], 'det_rout',[], 'det_side',side_det, ...
         'det_offset_mm',det_offset_sq_mm, ...
         'p_mode','sampled_only'), ...

  struct('name','7) Gaussian src + Gaussian-weighted det (overlap, cap-fix)', ...
         'type','gaussian_overlap', ...
         'w_mm',w_mm, 'R_det_mm',R_det_g, 'Nw_cap',Nw_cap, ...
         'p_mode','analytic_gaussian_cap') ...
};

% Apply CasesToRun mask
keepIdx = find(logical(CasesToRun(:).'));
Cases = CasesAll(keepIdx);
nC = numel(Cases);

if nC==0
    error('No cases selected. Set CasesToRun(i)=1 for at least one case.');
end

%% ========================================================================
%  4) Precompute geometry points (Fig 1) + p(rho) (Fig 2 + mixing)
% ========================================================================

geom_src_pts   = cell(nC,1);
geom_det_pts   = cell(nC,1);

rho_mid_store   = cell(nC,1);
rho_edges_store = cell(nC,1);
Pmix_store      = cell(nC,1);   % p(rho) used for mixing (normalized)
Ptheo_store     = cell(nC,1);   % theoretical overlay (if available)
Ptheo_label     = cell(nC,1);
rho_samp_store  = cell(nC,1);   % sampled rho for Fig 2 sanity curve

A_det_store     = nan(nC,1);    % detector area used in SIA (or effective area for Gaussian)

for c = 1:nC
    C = Cases{c};

    switch C.type
        case 'uniform_shapes'
            [Ps_plot, Pd_plot] = sample_geometry_points_uniform(C, max(Ngeom_plot,1000));
            geom_src_pts{c} = Ps_plot(1:min(Ngeom_plot,size(Ps_plot,1)),:);
            geom_det_pts{c} = Pd_plot(1:min(Ngeom_plot,size(Pd_plot,1)),:);

            rho_max = support_rhomax_uniform(C);
            rho_edges = linspace(0, rho_max, nbins_rho+1);
            rho_mid   = (rho_edges(1:end-1) + rho_edges(2:end))/2;
            drho      = rho_edges(2) - rho_edges(1);

            switch C.p_mode
                case 'analytic_circle'
                    Pmix = p_rho_circle(rho_mid, 2*C.det_rout); % df=2R
                    Ptheo = Pmix; Ptheo_label{c} = 'analytic';
                case 'analytic_core_total'
                    Pmix = p_rho_core_total_mid(rho_mid, C.src_rout, C.det_rout);
                    Ptheo = Pmix; Ptheo_label{c} = 'analytic';
                case 'analytic_DRS'
                    pfun = p_rho_DRS(norm(C.det_offset_mm), C.src_rout, C.det_rout);
                    Pmix = pfun(rho_mid);
                    Ptheo = Pmix; Ptheo_label{c} = 'semi-analytic (1D integral)';
                case 'sampled_only'
                    rho_hi = sample_rho_geometry_uniform(C, Ngeom_pairs_mix);
                    Pmix = hist_pdf_from_samples_mid(rho_hi, rho_mid);
                    Ptheo = []; Ptheo_label{c} = 'sampled-only';
                otherwise
                    error('Unknown p_mode for uniform_shapes.');
            end

            Pmix = sanitize_pdf(Pmix);
            Z = sum(Pmix) * drho;
            if Z > 0, Pmix = Pmix / Z; end

            rho_samp = sample_rho_geometry_uniform(C, Ngeom_pairs_prho);
            A_det_store(c) = detector_area_uniform(C.det_shape, C.det_rin, C.det_rout);

        case 'square_square'
            [Ps_plot, Pd_plot] = sample_geometry_points_square(C, max(Ngeom_plot,1000));
            geom_src_pts{c} = Ps_plot(1:min(Ngeom_plot,size(Ps_plot,1)),:);
            geom_det_pts{c} = Pd_plot(1:min(Ngeom_plot,size(Pd_plot,1)),:);

            rho_max = hypot(C.src_side + C.det_side, C.src_side + C.det_side);
            rho_edges = linspace(0, rho_max, nbins_rho+1);
            rho_mid   = (rho_edges(1:end-1) + rho_edges(2:end))/2;
            drho      = rho_edges(2) - rho_edges(1);

            rho_hi = sample_rho_geometry_square(C, Ngeom_pairs_mix);
            Pmix = hist_pdf_from_samples_mid(rho_hi, rho_mid);
            Pmix = sanitize_pdf(Pmix);
            Z = sum(Pmix) * drho;
            if Z > 0, Pmix = Pmix / Z; end

            Ptheo = []; Ptheo_label{c} = 'sampled-only';
            rho_samp = sample_rho_geometry_square(C, Ngeom_pairs_prho);

            A_det_store(c) = C.det_side^2;

        case 'gaussian_overlap'
            w  = C.w_mm;
            Rd = C.R_det_mm;
            Rsrc_cap = Rd + C.Nw_cap * w;
            rho_max  = Rsrc_cap + Rd;

            [Ps_plot, Pd_plot] = sample_geometry_points_gaussian_overlap(max(Ngeom_plot,1000), w, Rsrc_cap, Rd, [0 0]);
            geom_src_pts{c} = Ps_plot(1:min(Ngeom_plot,size(Ps_plot,1)),:);
            geom_det_pts{c} = Pd_plot(1:min(Ngeom_plot,size(Pd_plot,1)),:);

            rho_edges = linspace(0, rho_max, nbins_rho+1);
            rho_mid   = (rho_edges(1:end-1) + rho_edges(2:end))/2;
            drho      = rho_edges(2) - rho_edges(1);

            Pmix = peff_gaussian_overlap_integral(rho_mid, drho, w, Rsrc_cap, Rd);
            Ptheo = Pmix;
            Ptheo_label{c} = sprintf('analytic integral (cap: R_{src}=R_{det}+%dw)', C.Nw_cap);

            rho_samp = sample_rho_gauss_trunc_pairs(Ngeom_pairs_prho, w, Rsrc_cap, Rd);

            A_det_store(c) = (pi*w^2/2) * (1 - exp(-2*Rd^2/w^2));

        otherwise
            error('Unknown case type.');
    end

    rho_mid_store{c}   = rho_mid;
    rho_edges_store{c} = rho_edges;
    Pmix_store{c}      = Pmix;

    Ptheo_store{c} = sanitize_pdf(Ptheo);
    if ~isempty(Ptheo_store{c})
        Zt = sum(Ptheo_store{c}) * drho;
        if Zt > 0, Ptheo_store{c} = Ptheo_store{c} / Zt; end
    end

    rho_samp_store{c} = rho_samp;
end

% Global rho_max for pencil trimming (since we reuse the same pencil run for all cases)
rho_max_all = 0;
for c = 1:nC
    rho_max_all = max(rho_max_all, rho_edges_store{c}(end));
end

%% ========================================================================
%  5) Allocate outputs
% ========================================================================

R_full = nan(nC, reps_per_cfg);
R_sia  = nan(nC, reps_per_cfg);

meanL_full = nan(nC, reps_per_cfg);
meanL_sia  = nan(nC, reps_per_cfg);
W1_L       = nan(nC, reps_per_cfg);

pL_full_mean = zeros(nC, nLbins);
pL_sia_mean  = zeros(nC, nLbins);

%% ========================================================================
%  6) Main loop (OPTION A): reps outer, pencil once per rep, FSD per case per rep
% ========================================================================

% Build pencil config ONCE
cfg_pencil = cfg0;
cfg_pencil.srctype = 'pencil';
cfg_pencil = rmfield_safe(cfg_pencil, 'srcparam1');
cfg_pencil = rmfield_safe(cfg_pencil, 'srcparam2');

% Progress steps: per rep => 1 pencil + nC FSD runs
totalSteps = reps_per_cfg * (nC + 1);
step = 0;

wb = waitbar(0, 'Starting...', 'Name','MCX progress');
cleanupObj = onCleanup(@() close_waitbar_safe(wb));

for r = 1:reps_per_cfg
    % ==========================================================
    % (1) Pencil run ONCE for this rep (shared by all scenarios)
    % ==========================================================
    step = step + 1;
    if ishandle(wb)
        waitbar(step/totalSteps, wb, sprintf('Rep %d/%d | Pencil (shared)', r, reps_per_cfg));
    end

    [~, detP] = mcxlab(cfg_pencil);

    % Pre-extract and trim pencil detections to GLOBAL rho_max_all and NA
    p_has = isfield(detP,'p') && ~isempty(detP.p) && isfield(detP,'data') && ~isempty(detP.data);
    if ~p_has
        rhoP_all = zeros(0,1);
        L_P_all  = zeros(0,1);
        wP_all   = zeros(0,1);
        dxP_all  = zeros(0,1);
        dyP_all  = zeros(0,1);
        det_data_all = zeros(4,0);
    else
        [dataP, posP_vox] = mask_by_radius_and_NA_mm(detP, cfg_pencil.unitinmm, cfg_pencil.srcpos(1:2), theta_acc_max, rho_max_all);

        if isempty(posP_vox)
            rhoP_all = zeros(0,1);
            L_P_all  = zeros(0,1);
            wP_all   = zeros(0,1);
            dxP_all  = zeros(0,1);
            dyP_all  = zeros(0,1);
            det_data_all = zeros(4,0);
        else
            posP_mm = double(posP_vox) * cfg_pencil.unitinmm;
            dxP_all = posP_mm(:,1) - src_center_mm(1);
            dyP_all = posP_mm(:,2) - src_center_mm(2);

            rhoP_all = hypot(dxP_all, dyP_all);
            L_P_all  = double(dataP(3,:)).';
            wP_all   = exp(-mu_a .* L_P_all);

            keep = isfinite(rhoP_all) & isfinite(L_P_all) & isfinite(wP_all) & ...
                   (rhoP_all>=0) & (rhoP_all<=rho_max_all) & (L_P_all>=0);

            rhoP_all = rhoP_all(keep);
            L_P_all  = L_P_all(keep);
            wP_all   = wP_all(keep);
            dxP_all  = dxP_all(keep);
            dyP_all  = dyP_all(keep);

            det_data_all = zeros(4, numel(L_P_all));
            det_data_all(2,:) = L_P_all.';   % L
            det_data_all(3,:) = dxP_all.';   % dx
            det_data_all(4,:) = dyP_all.';   % dy
        end
    end

    % ==========================================================
    % (2) Per-case processing for this rep:
    %     - FSD: run MCX per case
    %     - SIA: reuse pencil detections (NO NEW MCX)
    % ==========================================================
    for c = 1:nC
        C = Cases{c};

        step = step + 1;
        if ishandle(wb)
            waitbar(step/totalSteps, wb, sprintf('Rep %d/%d | Case %d/%d FSD', r, reps_per_cfg, c, nC));
        end

        rho_mid   = rho_mid_store{c};
        rho_edges = rho_edges_store{c};
        Pbin      = Pmix_store{c};
        drho      = rho_edges(2) - rho_edges(1);

        annA = pi*(rho_edges(2:end).^2 - rho_edges(1:end-1).^2);
        annA = double(annA(:));
        mixR = (double(Pbin(:)) .* double(drho)) ./ max(double(annA(:)), eps);

        A_det_use = A_det_store(c);

        % ---------------- Build MCX config for THIS case (FSD) ----------------
        cfg_full = cfg0;
        det_center_mm = src_center_mm;
        if isfield(C,'det_offset_mm')
            det_center_mm = src_center_mm + double(C.det_offset_mm(:)).';
        end

        switch C.type
            case 'uniform_shapes'
                cfg_full.srctype = 'disk';
                cfg_full = rmfield_safe(cfg_full, 'srcparam1');
                cfg_full = rmfield_safe(cfg_full, 'srcparam2');

                if strcmp(C.src_shape,'disk')
                    cfg_full.srcparam1 = [C.src_rout, 0] / cfg_full.unitinmm;      % [R rin] vox
                elseif strcmp(C.src_shape,'annulus')
                    cfg_full.srcparam1 = [C.src_rout, C.src_rin] / cfg_full.unitinmm;
                else
                    error('Unsupported source shape for uniform_shapes.');
                end

            case 'square_square'
                cfg_full.srctype = 'planar';
                cfg_full = rmfield_safe(cfg_full, 'srcparam1');
                cfg_full = rmfield_safe(cfg_full, 'srcparam2');

                side_vox = C.src_side / cfg_full.unitinmm;
                corner_vox = [double(cfg0.srcpos(1)) - side_vox/2, ...
                              double(cfg0.srcpos(2)) - side_vox/2, ...
                              double(cfg0.srcpos(3))];
                cfg_full.srcpos    = corner_vox;
                cfg_full.srcparam1 = [side_vox 0 0 0];
                cfg_full.srcparam2 = [0 side_vox 0 0];
                cfg_full.srcdir    = [0 0 1];

            case 'gaussian_overlap'
                cfg_full.srctype   = 'gaussian';
                cfg_full.srcdir    = [0 0 1 0];
                cfg_full.srcparam1 = [C.w_mm / cfg_full.unitinmm 0 0 0];
                cfg_full.srcparam2 = [0 0 0 0];

            otherwise
                error('Unknown case type in MCX build.');
        end

        % ==========================================================
        % FSD: run + detector mask (+ optional detector weighting)
        % ==========================================================
        [~, detF] = mcxlab(cfg_full);

        pL_F = zeros(1, nLbins);

        switch C.type
            case 'uniform_shapes'
                [dataF, posF_vox] = mask_by_detector_general_mm(detF, cfg_full.unitinmm, det_center_mm, theta_acc_max, C.det_shape, C.det_rin, C.det_rout);
                if isempty(posF_vox)
                    R_full(c,r) = 0;
                else
                    L_F = double(dataF(3,:));
                    w_F = exp(-mu_a .* L_F);
                    R_full(c,r) = sum(w_F) / Nlaunch;
                    pL_F = weighted_pdf_1d_accum(L_F, w_F, edges_L);
                end

            case 'square_square'
                [dataF, posF_vox] = mask_by_square_and_NA_mm(detF, cfg_full.unitinmm, det_center_mm, C.det_side, theta_acc_max);
                if isempty(posF_vox)
                    R_full(c,r) = 0;
                else
                    L_F = double(dataF(3,:));
                    w_F = exp(-mu_a .* L_F);
                    R_full(c,r) = sum(w_F) / Nlaunch;
                    pL_F = weighted_pdf_1d_accum(L_F, w_F, edges_L);
                end

            case 'gaussian_overlap'
                if ~isfield(detF,'p') || isempty(detF.p)
                    R_full(c,r) = 0;
                else
                    pos_xy_vox = detF.p(:,1:2);
                    pos_xy_mm  = double(pos_xy_vox) * cfg_full.unitinmm;
                    dx = pos_xy_mm(:,1) - src_center_mm(1);
                    dy = pos_xy_mm(:,2) - src_center_mm(2);
                    r_det = hypot(dx, dy);

                    ang = acos(-double(detF.data(9,:)));
                    mask = (ang.' <= theta_acc_max) & (r_det <= C.R_det_mm);

                    if ~any(mask)
                        R_full(c,r) = 0;
                    else
                        L_F  = double(detF.data(3, mask));
                        wAbs = exp(-mu_a .* L_F);
                        wDet = exp(-2*(r_det(mask).^2)/(C.w_mm^2));
                        wFin = wAbs(:) .* wDet(:);

                        R_full(c,r) = sum(wFin) / Nlaunch;
                        pL_F = weighted_pdf_1d_accum(L_F(:), wFin(:), edges_L);
                    end
                end

            otherwise
                error('Unknown case type in FSD.');
        end

        % ==========================================================
        % SIA: reuse shared pencil detections from THIS rep (NO MCX call)
        % ==========================================================
        pL_S = zeros(1, nLbins);

        if isempty(rhoP_all)
            R_sia(c,r) = 0;
        else
            % ---- SIA reflectance (absolute) ----
            binR = discretize(rhoP_all, rho_edges);
            okR  = ~isnan(binR);
            if ~any(okR)
                R_sia(c,r) = 0;
            else
                sumsR = accumarray(double(binR(okR)), double(wP_all(okR)), [numel(rho_mid), 1], @sum, 0);
                R_sia(c,r) = (A_det_use / Nlaunch) * sum(double(sumsR(:)) .* double(mixR(:)));
            end

            % ---- SIA p(L) (shape-only) ----
            % Uses same det_data_all; internal trimming uses this case's rho_edges
            pL_S = weighted_pdf_L_sia_accum(det_data_all, mu_a, Nlaunch, A_det_use, rho_edges, edges_L, Pbin, drho);
        end

        % ==========================================================
        % Metrics (per rep)
        % ==========================================================
        meanL_full(c,r) = sum(L_mid .* pL_F) * dL;
        meanL_sia (c,r) = sum(L_mid .* pL_S) * dL;

        cdfF = cumsum(pL_F) * dL;
        cdfS = cumsum(pL_S) * dL;
        W1_L(c,r) = sum(abs(cdfF - cdfS)) * dL;

        % Running mean PDFs (per case)
        pL_full_mean(c,:) = pL_full_mean(c,:) + (pL_F - pL_full_mean(c,:)) / r;
        pL_sia_mean (c,:) = pL_sia_mean (c,:) + (pL_S - pL_sia_mean (c,:)) / r;

        fprintf('rep %2d/%2d | %s | R_full=%.6g | R_sia=%.6g | <L>F=%.3f | <L>S=%.3f | W1=%.4g\n', ...
            r, reps_per_cfg, C.name, R_full(c,r), R_sia(c,r), meanL_full(c,r), meanL_sia(c,r), W1_L(c,r));
    end
end

if ishandle(wb), waitbar(1, wb, 'Done'); end

%% ========================================================================
%  7) Summary (per configuration)
% ========================================================================

fprintf('\n====================== SUMMARY ======================\n');
for c = 1:nC
    mRF = mean(R_full(c,:), 'omitnan'); sRF = std(R_full(c,:), 0, 'omitnan');
    mRS = mean(R_sia (c,:), 'omitnan'); sRS = std(R_sia (c,:), 0, 'omitnan');

    mLF = mean(meanL_full(c,:), 'omitnan'); sLF = std(meanL_full(c,:), 0, 'omitnan');
    mLS = mean(meanL_sia (c,:), 'omitnan'); sLS = std(meanL_sia (c,:), 0, 'omitnan');

    mW1 = mean(W1_L(c,:), 'omitnan');
    mW1_med = median(W1_L(c,:), 'omitnan');

    fprintf('\n%s\n', Cases{c}.name);
    fprintf('  Reflectance: Full = %.6g ± %.3g | SIA = %.6g ± %.3g | Δmean = %.3g\n', ...
        mRF, sRF, mRS, sRS, (mRS - mRF));
    fprintf('  Mean pathlength <L> [mm] (shape-only p(L)): Full = %.6g ± %.3g | SIA = %.6g ± %.3g | Δmean = %.3g\n', ...
        mLF, sLF, mLS, sLS, (mLS - mLF));
    fprintf('  W1(CDF) [mm]: mean = %.4g | median = %.4g\n', mW1, mW1_med);

    if strcmp(Cases{c}.type,'gaussian_overlap')
        w  = Cases{c}.w_mm;
        Rd = Cases{c}.R_det_mm;
        Rsrc = Rd + Cases{c}.Nw_cap*w;
        fprintf('  Gaussian details: w=%.3f mm, R_det=%.3f mm, R_src_cap=%.3f mm, A_det_eff=%.6g mm^2\n', ...
            w, Rd, Rsrc, A_det_store(c));
    end
end

%% ========================================================================
%  8) FIG 1 — Geometry: sampled input/output points (source-centered)
% ========================================================================

figure('Color','w');
t1 = tiledlayout(nC,1,'TileSpacing','compact','Padding','compact');

for c = 1:nC
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on'); axis(ax,'equal');

    Ps = geom_src_pts{c};
    Pd = geom_det_pts{c};

    Ns = min(500, size(Ps,1));
    Nd = min(500, size(Pd,1));

    scatter(ax, Ps(1:Ns,1), Ps(1:Ns,2), 6, 'filled', 'DisplayName','Input (source)');
    scatter(ax, Pd(1:Nd,1), Pd(1:Nd,2), 6, 'filled', 'DisplayName','Output (detector)');

    xlabel(ax,'x [mm]');
    ylabel(ax,'y [mm]');
    title(ax, Cases{c}.name, 'Interpreter','none');

    Rmax = max([hypot(Ps(:,1),Ps(:,2)); hypot(Pd(:,1),Pd(:,2)); 1e-3]);
    lim  = 1.10 * Rmax;
    xlim(ax, [-lim lim]);
    ylim(ax, [-lim lim]);

    if c == 1
        legend(ax,'Location','best');
    end
end

%% ========================================================================
%  9) FIG 2 — p(rho): sampled sanity vs theoretical (when available)
% ========================================================================

figure('Color','w');
t2 = tiledlayout(nC,1,'TileSpacing','compact','Padding','compact');

for c = 1:nC
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    rho_mid  = rho_mid_store{c};
    rho_samp = rho_samp_store{c};

    edges = linspace(0, max(rho_mid), Nbins_prho_plot+1);
    counts = histcounts(rho_samp, edges);
    binw = diff(edges);
    pdf_samp = counts ./ (sum(counts).*binw + eps);
    rho_cent = (edges(1:end-1)+edges(2:end))/2;

    plot(ax, rho_cent, pdf_samp, '-', 'LineWidth', 2.2, 'DisplayName','Sampled');

    Ptheo = Ptheo_store{c};
    if ~isempty(Ptheo)
        plot(ax, rho_mid, Ptheo, '--', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('Theoretical (%s)', Ptheo_label{c}));
    end

    xlabel(ax,'\rho [mm]');
    ylabel(ax,'p(\rho)');
    title(ax, Cases{c}.name, 'Interpreter','none');

    if c == 1
        legend(ax,'Location','best');
    end
end

%% ========================================================================
% 10) FIG 3 — Reflectance distributions (hist overlays)
% ========================================================================

figure('Color','w');
t3 = tiledlayout(nC,1,'TileSpacing','compact','Padding','compact');

for c = 1:nC
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    rf = R_full(c,:); rf = rf(isfinite(rf));
    rs = R_sia (c,:); rs = rs(isfinite(rs));

    xmin = min([rf(:); rs(:)]);
    xmax = max([rf(:); rs(:)]);
    if ~isfinite(xmin) || ~isfinite(xmax) || xmin==xmax
        xmin = 0; xmax = 1;
    end

    edgesR = linspace(max(0,xmin), xmax, 40);

    histogram(ax, rf, edgesR, 'Normalization','pdf', 'FaceAlpha',0.35, 'EdgeColor','none', 'DisplayName','Full S-D');
    histogram(ax, rs, edgesR, 'Normalization','pdf', 'FaceAlpha',0.35, 'EdgeColor','none', 'DisplayName','SIA');

    xlabel(ax,'Reflectance R');
    ylabel(ax,'PDF');
    title(ax, Cases{c}.name, 'Interpreter','none');

    if c == 1
        legend(ax,'Location','best');
    end
end

%% ========================================================================
% 11) FIG 4 — Mean detected pathlength PDFs p(L) (avg over reps)
% ========================================================================

figure('Color','w');
t4 = tiledlayout(nC,1,'TileSpacing','compact','Padding','compact');

for c = 1:nC
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    plot(ax, L_mid, pL_full_mean(c,:), '-',  'LineWidth', 2.0, 'DisplayName','Full S-D');
    plot(ax, L_mid, pL_sia_mean (c,:), '--', 'LineWidth', 2.0, 'DisplayName','SIA');

    xlabel(ax,'Pathlength L [mm]');
    ylabel(ax,'p(L) (absorption-weighted, normalized)');
    title(ax, Cases{c}.name, 'Interpreter','none');

    if c == 1
        legend(ax,'Location','best');
    end
end

%% ========================================================================
%  Local functions
% ========================================================================

function rho_max = support_rhomax_uniform(C)
    if strcmp(C.src_shape,'disk') && strcmp(C.det_shape,'disk') && all(C.det_offset_mm==0) && (C.src_rout==C.det_rout)
        rho_max = 2*C.det_rout;
        return;
    end
    if all(C.det_offset_mm==0)
        rho_max = C.src_rout + C.det_rout;
    else
        rho_max = norm(C.det_offset_mm) + (C.src_rout + C.det_rout);
    end
end

function A = detector_area_uniform(shape, rin, rout)
    switch shape
        case 'disk'
            A = pi*rout^2;
        case 'annulus'
            A = pi*(rout^2 - rin^2);
        otherwise
            error('Unknown uniform detector shape');
    end
end

function [data_cut, pos_cut_vox] = mask_by_radius_and_NA_mm(det, unitinmm, srcpos_vox_xy, theta_acc_max, r_max_mm)
    if ~isfield(det,'p') || isempty(det.p)
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end
    pos_xy_vox = det.p(:,1:2);
    d_vox = double(pos_xy_vox) - double(srcpos_vox_xy);
    r_mm = hypot(d_vox(:,1), d_vox(:,2)) * unitinmm;

    ang = acos(-double(det.data(9,:)));
    mask = (r_mm <= r_max_mm) & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask,:);
    else
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function [data_cut, pos_cut_vox] = mask_by_detector_general_mm(det, unitinmm, det_center_mm, theta_acc_max, det_shape, det_rin, det_rout)
    if ~isfield(det,'p') || isempty(det.p)
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end
    pos_xy_vox = det.p(:,1:2);
    pos_xy_mm  = double(pos_xy_vox) * unitinmm;
    r_det_mm   = hypot(pos_xy_mm(:,1)-det_center_mm(1), pos_xy_mm(:,2)-det_center_mm(2));

    ang = acos(-double(det.data(9,:)));

    switch det_shape
        case 'disk'
            radial_ok = (r_det_mm <= det_rout);
        case 'annulus'
            radial_ok = (r_det_mm >= det_rin) & (r_det_mm <= det_rout);
        otherwise
            error('Unsupported det_shape in mask_by_detector_general_mm.');
    end

    mask = radial_ok & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask,:);
    else
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function [data_cut, pos_cut_vox] = mask_by_square_and_NA_mm(det, unitinmm, det_center_mm, det_side_mm, theta_acc_max)
    if ~isfield(det,'p') || isempty(det.p)
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
        return;
    end
    pos_xy_vox = det.p(:,1:2);
    pos_xy_mm  = double(pos_xy_vox) * unitinmm;

    dx = pos_xy_mm(:,1) - det_center_mm(1);
    dy = pos_xy_mm(:,2) - det_center_mm(2);
    half = det_side_mm/2;

    in_square = (abs(dx) <= half) & (abs(dy) <= half);

    ang = acos(-double(det.data(9,:)));
    mask = in_square & (ang.' <= theta_acc_max);

    if any(mask)
        data_cut = det.data(:, mask);
        pos_cut_vox = pos_xy_vox(mask,:);
    else
        data_cut = det.data(:, false(1,0));
        pos_cut_vox = zeros(0,2);
    end
end

function pdf = weighted_pdf_1d_accum(x, w, edges)
    nb = numel(edges) - 1;
    pdf = zeros(1, nb);
    if isempty(x), return; end

    x = double(x(:));
    w = double(w(:));

    idx = discretize(x, edges);
    ok  = ~isnan(idx);
    if ~any(ok), return; end

    idx = double(idx(ok));
    w   = double(w(ok));

    wk = accumarray(idx(:), w(:), [nb, 1], @sum, 0);

    dx = edges(2) - edges(1);
    Z  = sum(wk) * dx;
    if Z > 0
        pdf = (wk / Z).';
    end
end

function pL = weighted_pdf_L_sia_accum(det_data, mu_a, launched_photons, A_det_use, edges_rho, edges_L, Pbin, drho)
    nR = numel(edges_rho) - 1;
    nL = numel(edges_L)   - 1;
    pL = zeros(1, nL);

    if isempty(det_data) || size(det_data,2)==0
        return;
    end

    rho = double(sqrt(det_data(3,:).^2 + det_data(4,:).^2));
    L   = double(det_data(2,:));
    w   = exp(-mu_a .* L);

    keep = (rho >= edges_rho(1)) & (rho <= edges_rho(end)) & ...
           (L   >= edges_L(1))   & (L   <= edges_L(end));
    if ~any(keep), return; end

    rho = rho(keep);
    L   = L(keep);
    w   = double(w(keep));

    ir = discretize(rho, edges_rho);
    iL = discretize(L,   edges_L);

    ok = ~isnan(ir) & ~isnan(iL);
    if ~any(ok), return; end

    ir = double(ir(ok));
    iL = double(iL(ok));
    w  = double(w(ok));

    W2 = accumarray([ir(:), iL(:)], w(:), [nR, nL], @sum, 0);

    ann = pi * (edges_rho(2:end).^2 - edges_rho(1:end-1).^2);
    ann = double(ann(:));

    R2 = W2 ./ (ann * launched_photons);

    PL = A_det_use * ((double(Pbin(:)) .* drho).' * R2);

    dL = edges_L(2) - edges_L(1);
    Z  = sum(PL) * dL;
    if Z > 0
        pL = PL / Z;
    end
end

function p = p_rho_circle(rho, df)
    rho = rho(:);
    x = rho./df;
    x = max(0, min(1, x));
    p = (16.*rho)./(pi*df.^2).*acos(x) - (16.*rho.^2)./(pi*df.^3).*sqrt(1 - x.^2);
    p(~isfinite(p)) = 0;
    p(rho<0 | rho>df) = 0;
end

function p = p_rho_core_total_mid(rho, r_in, R)
    rho = rho(:);
    p = zeros(size(rho));
    for i = 1:numel(rho)
        s = rho(i);
        if s < 0 || s > (R + r_in)
            p(i) = 0;
        elseif s <= (R - r_in)
            p(i) = (2*s)/(R^2);
        else
            a1 = (r_in^2 + s^2 - R^2)/(2*r_in*s);
            a2 = (R^2 + s^2 - r_in^2)/(2*R*s);
            alpha = acos(max(-1,min(1,a1)));
            beta  = acos(max(-1,min(1,a2)));
            p(i) = (s/pi)*((2*alpha - sin(2*alpha))/R^2 + (2*beta - sin(2*beta))/r_in^2);
        end
    end
    p(~isfinite(p)) = 0;
    p(p<0) = 0;
end

function P = sanitize_pdf(P)
    if isempty(P), return; end
    P = double(P(:));
    P(~isfinite(P)) = 0;
    P(P<0) = 0;
end

function [Ps, Pd] = sample_geometry_points_uniform(C, N)
    Ps = sample_region_uniform(C.src_shape, N, C.src_rin, C.src_rout);
    Pd = sample_region_uniform(C.det_shape, N, C.det_rin, C.det_rout) + double(C.det_offset_mm(:)).';
end

function rho = sample_rho_geometry_uniform(C, Npairs)
    Ps = sample_region_uniform(C.src_shape, Npairs, C.src_rin, C.src_rout);
    Pd = sample_region_uniform(C.det_shape, Npairs, C.det_rin, C.det_rout) + double(C.det_offset_mm(:)).';
    rho = hypot(Ps(:,1) - Pd(:,1), Ps(:,2) - Pd(:,2));
end

function P = sample_region_uniform(shape, N, rin, rout)
    switch shape
        case 'disk'
            th = 2*pi*rand(N,1);
            r  = rout*sqrt(rand(N,1));
            P  = [r.*cos(th), r.*sin(th)];
        case 'annulus'
            th = 2*pi*rand(N,1);
            r  = sqrt((rout^2 - rin^2)*rand(N,1) + rin^2);
            P  = [r.*cos(th), r.*sin(th)];
        otherwise
            error('Unknown uniform region shape.');
    end
end

function [Ps, Pd] = sample_geometry_points_square(C, N)
    Ps = sample_square(N, C.src_side);
    Pd = sample_square(N, C.det_side) + double(C.det_offset_mm(:)).';
end

function rho = sample_rho_geometry_square(C, Npairs)
    Ps = sample_square(Npairs, C.src_side);
    Pd = sample_square(Npairs, C.det_side) + double(C.det_offset_mm(:)).';
    rho = hypot(Ps(:,1) - Pd(:,1), Ps(:,2) - Pd(:,2));
end

function P = sample_square(N, side_mm)
    half = side_mm/2;
    P = (rand(N,2) - 0.5) * (2*half);
end

function P_mid = hist_pdf_from_samples_mid(rho_samp, rho_mid)
    rho_samp = rho_samp(:);
    rho_mid  = rho_mid(:);
    dr = rho_mid(2) - rho_mid(1);
    edges = [0; rho_mid + dr/2];
    counts = histcounts(rho_samp, edges);
    counts = double(counts(:));
    P_mid = counts ./ (sum(counts)*dr + eps);
end

function p_rho = p_rho_DRS(r_sep, a, b)
    Ntheta = 2048;
    theta  = linspace(0, 2*pi, Ntheta);
    cth    = cos(theta);
    p_rho  = @(rho) pdf_eval_DRS(rho, r_sep, a, b, theta, cth);
end

function f = pdf_eval_DRS(rho, r, a, b, theta, cth)
    rho = rho(:);
    f   = zeros(size(rho));

    dmin  = r - (a + b);
    dmax  = r + (a + b);
    normC = (pi^2) * (a^2) * (b^2);

    for i = 1:numel(rho)
        d = rho(i);
        if (d < dmin) || (d > dmax) || (d < 0)
            f(i) = 0;
            continue;
        end

        s  = sqrt(r^2 + d^2 + 2*r*d*cth);
        Ac = overlap_area(s, a, b);
        I  = trapz(theta, Ac);

        f(i) = (d / normC) * I;
    end
end

function A = overlap_area(s, a, b)
    s = max(s, 0);
    A = zeros(size(s));

    noOverlap = (s >= (a + b));
    contain   = (s <= abs(a - b));
    partial   = ~(noOverlap | contain);

    A(contain) = pi * (min(a,b)^2);

    sp = s(partial);
    c1 = (sp.^2 + a^2 - b^2) ./ (2*sp*a);
    c2 = (sp.^2 + b^2 - a^2) ./ (2*sp*b);
    c1 = min(max(c1, -1), 1);
    c2 = min(max(c2, -1), 1);

    term = (-sp + a + b) .* (sp + a - b) .* (sp - a + b) .* (sp + a + b);
    term = max(term, 0);

    A(partial) = a^2 .* acos(c1) + b^2 .* acos(c2) - 0.5 .* sqrt(term);
end

function [Ps, Pd] = sample_geometry_points_gaussian_overlap(N, w_mm, R_s, R_d, det_offset_mm)
    inv_r = @(u, Rtr) (w_mm/sqrt(2)) * sqrt( -log( 1 - u*(1 - exp(-2*Rtr^2/w_mm^2)) ) );

    u1  = rand(N,1);  th1 = 2*pi*rand(N,1);
    u2  = rand(N,1);  th2 = 2*pi*rand(N,1);

    rs  = inv_r(u1, R_s);
    rd  = inv_r(u2, R_d);

    Ps  = [rs.*cos(th1), rs.*sin(th1)];
    Pd  = [rd.*cos(th2), rd.*sin(th2)] + det_offset_mm;
end

function rho = sample_rho_gauss_trunc_pairs(Npairs, w_mm, R_s, R_d)
    inv_r = @(u, Rtr) (w_mm/sqrt(2)) * sqrt( -log( 1 - u*(1 - exp(-2*Rtr^2/w_mm^2)) ) );

    u1  = rand(Npairs,1);  th1 = 2*pi*rand(Npairs,1);
    u2  = rand(Npairs,1);  th2 = 2*pi*rand(Npairs,1);

    rs  = inv_r(u1, R_s);
    rd  = inv_r(u2, R_d);

    xs  = rs .* cos(th1);  ys = rs .* sin(th1);
    xd  = rd .* cos(th2);  yd = rd .* sin(th2);

    rho = hypot(xs - xd, ys - yd);
end

function p = peff_gaussian_overlap_integral(rho_mid, drho, w_mm, R_s, R_d)
    rho_mid = rho_mid(:);

    Z = @(R) (pi*w_mm^2/2) * (1 - exp(-2*R.^2/w_mm^2));
    Zs = Z(R_s);  Zd = Z(R_d);

    qs = @(r) exp(-2*(r.^2)/(w_mm^2)) / Zs;
    qd = @(r) exp(-2*(r.^2)/(w_mm^2)) / Zd;

    Nk   = 3000;
    kmax = 350;         % [1/mm]
    k    = linspace(0, kmax, Nk);

    Nr = 4500;
    rS = linspace(0, R_s, Nr);
    rD = linspace(0, R_d, Nr);

    qs_r = qs(rS) .* rS;
    qd_r = qd(rD) .* rD;

    Qs = zeros(size(k));
    Qd = zeros(size(k));

    for i = 1:numel(k)
        ki = k(i);
        Qs(i) = 2*pi * trapz(rS, qs_r .* besselj(0, ki*rS));
        Qd(i) = 2*pi * trapz(rD, qd_r .* besselj(0, ki*rD));
    end

    p = zeros(numel(rho_mid),1);
    for j = 1:numel(rho_mid)
        rho = rho_mid(j);
        integrand = (Qs .* Qd) .* besselj(0, k*rho) .* k;
        p(j) = rho * trapz(k, integrand);
    end

    p(~isfinite(p)) = 0;
    p(p<0) = 0;
    p = p ./ (sum(p)*drho + eps);
end

function s = rmfield_safe(s, f)
    if isfield(s,f), s = rmfield(s,f); end
end

function close_waitbar_safe(wb)
    if ~isempty(wb) && ishandle(wb)
        close(wb);
    end
end