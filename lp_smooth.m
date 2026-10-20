function [S, info] = lp_smooth(X, varargin)
% LP_SMOOTH  Log-penalty-based multiscale image smoothing.
%
%   Reference implementation of:
%     "Log-Penalty-Based Multiscale Image Smoothing for Edge Preservation
%      and Texture Suppression"
%
%   The method performs edge-preserving smoothing: it removes textures and
%   noise while keeping structural edges sharp, by minimizing a non-convex
%   energy with a logarithmic penalty and a multi-step gradient operator,
%   solved by an FFT-accelerated ADMM.
%
%   Usage:
%     S = lp_smooth(image)                    % one image in, one image out
%     S = lp_smooth(image, 'lambda',0.02)     % override any parameter
%     [S, info] = lp_smooth(image, ...)       % info = convergence logs
%
%   Input:
%     image   path to an image file, OR an image matrix (uint8 / double in
%             [0,1], grayscale or RGB).
%   Output:
%     S       smoothed image, double in [0,1], same size and channels.
%     info    struct: info.RE (relative error per iteration), info.objective.
%
%   Parameters (defaults are the optima reported in the paper):
%     'lambda'   smoothing strength.   default 0.01,  suggested [0.005, 0.03]
%     'kappa'    gradient sensitivity. default 0.30,  suggested [0.15, 0.40]
%     'gamma'    texture gain.         default 3.0,   suggested [2.0, 5.0]
%     'rho'      ADMM penalty.         default 0.5,   suggested [0.2, 0.7]
%     'max_iter' ADMM iterations.      default 20,    suggested [10, 30]
%
%   The four parameters above are the only tunable knobs; all other design
%   choices (directional grouping, softmax weighting, logarithmic penalty,
%   step sizes {1,2,4}) are fixed to the configuration evaluated in the paper.
%
%   Dependencies: MATLAB + Image Processing Toolbox.
%   Authors: Jianwu Long, Hao Luo, Shaoyi Wang (CQUT).

    % ---------- parameter parsing ----------
    p = inputParser;
    addParameter(p, 'lambda',   0.01,  @(x) isscalar(x) && x>0);
    addParameter(p, 'kappa',    0.30,  @(x) isscalar(x) && x>0);
    addParameter(p, 'gamma',    3.0,   @(x) isscalar(x) && x>=0);
    addParameter(p, 'rho',      0.5,   @(x) isscalar(x) && x>0);
    addParameter(p, 'max_iter', 20,    @(x) isscalar(x) && x>=1);
    parse(p, varargin{:});
    o = p.Results;

    % ---------- load image if a path was given ----------
    if ischar(X) || isstring(X)
        X = im2double(imread(char(X)));
    elseif ~isa(X, 'double')
        X = im2double(X);
    end

    scales    = [1, 2, 4];
    num_scales = numel(scales);
    [rows, cols, cha] = size(X);
    sizeI2D = [rows, cols];

    % ---------- 1. directional-grouped wavelet texture prior ----------
    TF = make_tight_frame_filters();
    R_prior = compute_texture_response(X, TF);

    % ---------- 2. ADMM initialization ----------
    S = X;
    R_aligned = zeros(rows, cols);
    dx = cell(1,num_scales); dy = cell(1,num_scales);
    zx = cell(1,num_scales); zy = cell(1,num_scales);
    Denormin = ones(rows, cols);

    for s = 1:num_scales
        k = scales(s);
        dx{s} = zeros(rows,cols,cha); dy{s} = zeros(rows,cols,cha);
        zx{s} = zeros(rows,cols,cha); zy{s} = zeros(rows,cols,cha);
        kx = zeros(1,k+1); kx(1)=1; kx(end)=-1;
        ky = zeros(k+1,1); ky(1)=1; ky(end)=-1;
        otfDx = psf2otf(kx, sizeI2D);
        otfDy = psf2otf(ky, sizeI2D);
        Denormin = Denormin + o.rho*(abs(otfDx).^2 + abs(otfDy).^2);
    end
    if cha > 1, Denormin = repmat(Denormin,[1,1,cha]); end
    fftX = fft2(X);

    ComputeGrad = @(im,k,d) im - circshift(im,-k,d);   % I(x) - I(x+k)
    ComputeDiv  = @(v,k,d)  v  - circshift(v, k,d);    % p(x) - p(x-k)

    RE_log = zeros(o.max_iter, 1);
    E_log  = zeros(o.max_iter, 1);

    % ---------- 3. ADMM main loop ----------
    for iter = 1:o.max_iter
        S_prev = S;

        % --- S-update (FFT) ---
        rhs_sum = zeros(rows, cols, cha);
        for s = 1:num_scales
            k = scales(s);
            rhs_sum = rhs_sum + ComputeDiv(dx{s}-zx{s},k,2) ...
                              + ComputeDiv(dy{s}-zy{s},k,1);
        end
        fftS = (fftX + fft2(o.rho*rhs_sum)) ./ Denormin;
        S = real(ifft2(fftS));

        % --- spatial weights (from the current iterate) ---
        if iter == 1
            R_aligned = imguidedfilter(R_prior, rgb2gray_custom(S), ...
                'NeighborhoodSize',[3 3],'DegreeOfSmoothing',1e-3);
        end

        % --- b-update (log-penalty prox) + z-update ---
        for s = 1:num_scales
            k = scales(s);
            gx = ComputeGrad(S,k,2);
            gy = ComputeGrad(S,k,1);
            gm = sqrt(sum(gx.^2,3) + sum(gy.^2,3));
            W = exp(-gm/o.kappa) .* (1 + o.gamma*R_aligned);
            alpha = (o.lambda/o.rho) * W;
            if cha > 1, alpha = repmat(alpha,[1,1,cha]); end

            dx{s} = prox_log(gx+zx{s}, alpha);
            dy{s} = prox_log(gy+zy{s}, alpha);
            zx{s} = zx{s} + gx - dx{s};
            zy{s} = zy{s} + gy - dy{s};
        end

        % --- convergence logging ---
        RE_log(iter) = norm(S(:)-S_prev(:)) / (norm(S(:))+eps);
        E_log(iter)  = compute_objective(X, S, o.lambda, o.gamma, R_aligned, o.kappa);
    end

    info = struct('RE', RE_log, 'objective', E_log);
end

% ===================== proximal operator of log(1+|t|) =====================

function u = prox_log(v, alpha)
% Closed-form proximal operator of mu*log(1+|t|).
    c = abs(v);
    disc = ((c+1)/2).^2 - alpha;
    root = sqrt(max(0, disc));
    u = sign(v) .* max(0, (c-1)/2 + root);
end

% ===================== objective value =====================

function E = compute_objective(X, S, lambda, gamma, R_aligned, kappa)
    scales = [1, 2, 4];
    fidelity = 0.5 * norm(S(:)-X(:))^2;
    reg = 0;
    for k = scales
        gx = S - circshift(S,-k,2);
        gy = S - circshift(S,-k,1);
        gm = sqrt(sum(gx.^2,3) + sum(gy.^2,3));
        W = exp(-gm/kappa) .* (1 + gamma*R_aligned);
        reg = reg + sum(W(:) .* log(gm(:) + 1));
    end
    E = fidelity + lambda * reg;
end

% ===================== tight-frame filter bank =====================

function filters = make_tight_frame_filters()
% 3x3 tight-frame (9-filter) bank. Index map: 1 = lowpass;
% 2,3 = horizontal group; 4,7 = vertical group; 5,6,8,9 = diagonal group.
    a1 = [1,2,1]/4; a2 = [sqrt(2),0,-sqrt(2)]/4; a3 = [-1,2,-1]/4;
    vecs = {a1,a2,a3};
    filters = cell(1,9);
    idx = 1;
    for i = 1:3
        for j = 1:3
            filters{idx} = vecs{i}' * vecs{j};
            idx = idx + 1;
        end
    end
end

% ===================== directional-grouped texture response =====================

function R = compute_texture_response(I, filters)
% Multi-scale directional RWDR response (Eqs. 9-14), fused by pixel-wise max.
    if size(I,3)>1, Ig = rgb2gray_custom(I); else, Ig = I; end
    [rows, cols] = size(Ig);
    levels = 3;
    R_acc = zeros(rows, cols);
    for l = 1:levels
        sc = 1/(2^(l-1));
        if l==1, cur = Ig; else, cur = imresize(Ig, sc, 'bicubic'); end
        Rs = compute_single_scale(cur, filters);
        if l==1, Rr = Rs; else, Rr = imresize(Rs,[rows,cols],'bicubic'); end
        R_acc = max(R_acc, Rr);
    end
    R = R_acc;
end

function R = compute_single_scale(img, filters)
% Directional-grouped RWDR at one scale with softmax weighting.
    sigma = 2.0; eps = 1e-5;
    [rows, cols] = size(img);
    nf = length(filters);
    C = zeros(rows, cols, nf);
    for k = 1:nf
        C(:,:,k) = imfilter(img, filters{k}, 'symmetric');
    end
    idx_H = [2,3]; idx_V = [4,7]; idx_D = [5,6,8,9];
    EH = imgaussfilt(sum(abs(C(:,:,idx_H)),3), sigma);
    EV = imgaussfilt(sum(abs(C(:,:,idx_V)),3), sigma);
    ED = imgaussfilt(sum(abs(C(:,:,idx_D)),3), sigma);
    WH = abs(imgaussfilt(sum(C(:,:,idx_H),3), sigma));
    WV = abs(imgaussfilt(sum(C(:,:,idx_V),3), sigma));
    WD = abs(imgaussfilt(sum(C(:,:,idx_D),3), sigma));

    % structural confidence alpha: 5-95 percentile-normalized variance
    ge = cat(3, EH, EV, ED);
    am = var(ge, 1, 3);
    amin = prctile(am(:), 5); amax = prctile(am(:), 95);
    an = (am - amin)/(amax - amin + 1e-6);
    an = max(0, min(1, an));

    % softmax over normalized energies, adaptive temperature T = exp(-alpha)
    se = EH + EV + ED + eps;
    T = exp(-an);
    eH = exp((EH./se)./T); eV = exp((EV./se)./T); eD = exp((ED./se)./T);
    es = eH + eV + eD + eps;
    wH = eH./es; wV = eV./es; wD = eD./es;

    R = (wH.*EH + wV.*EV + wD.*ED) ./ (wH.*WH + wV.*WV + wD.*WD + eps);

    % normalize to [0,1] (99th-percentile robust scaling)
    sp = R(1:10:end);
    mx = prctile(sp(:), 99); mn = min(sp(:));
    if mx <= mn, mx = mn + 1e-6; end
    R = max(mn, min(R, mx));
    R = (R - mn)/(mx - mn);
    R = imgaussfilt(R, 0.8);
end

% ===================== utilities =====================

function g = rgb2gray_custom(I)
    if size(I,3) == 3
        g = 0.2989*I(:,:,1) + 0.5870*I(:,:,2) + 0.1140*I(:,:,3);
    else
        g = I;
    end
end
