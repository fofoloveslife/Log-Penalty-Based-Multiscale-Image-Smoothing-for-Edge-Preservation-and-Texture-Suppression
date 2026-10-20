%% RUN_SMOOTHING  Smooth the image in this folder: just run it.
%
%   Put an image named "image" (e.g. image.png / image.jpg / image.bmp)
%   in the SAME folder as this script, then run:
%
%     run_smoothing
%
%   It smooths that image with the paper-optimal parameters, shows the
%   input and the smoothed result side by side, and saves the output as
%   'smoothed_image.png' in this folder.
%
%   Optional overrides:
%     run_smoothing('other_image.png')                 % smooth a specific file
%     run_smoothing([], 'lambda', 0.02)                % custom parameters
%     S = run_smoothing;                               % also return the result
%
%   Parameters (optional, key-value):
%     'lambda'   smoothing strength,   default 0.01,  range [0.005, 0.03]
%     'kappa'    gradient sensitivity, default 0.30,  range [0.15, 0.40]
%     'gamma'    texture gain,         default 3.0,   range [2.0, 5.0]
%     'rho'      ADMM penalty,         default 0.5,   range [0.2, 0.7]
%     'max_iter' ADMM iterations,      default 20,    range [10, 30]

function S = run_smoothing(image_path, varargin)

if nargin < 1 || isempty(image_path)
    image_path = find_local_image();
end

% ---- smooth (one call) ----
fprintf('Smoothing %s ...\n', image_path);
tic;
S = lp_smooth(image_path, varargin{:});
fprintf('Done in %.1f s.\n', toc);

% ---- display input vs output ----
img = im2double(imread(image_path));
figure('Name', 'Smoothing Result', 'Color', 'w', 'Position', [80, 80, 900, 400]);
subplot(1, 2, 1); imshow(img); title('Input', 'FontSize', 13);
subplot(1, 2, 2); imshow(S);   title('Smoothed', 'FontSize', 13);

% ---- save output ----
out_name = 'smoothed_image.png';
imwrite(S, out_name);
fprintf('Saved: %s\n', fullfile(pwd, out_name));

end

function p = find_local_image()
% Find an image named "image.*" in this folder; fall back to the first
% image file found; otherwise synthesize a test image.
here = fileparts(mfilename('fullpath'));
exts = {'png', 'jpg', 'jpeg', 'bmp', 'tif', 'tiff'};

% 1) look for "image.<ext>"
for i = 1:numel(exts)
    cand = fullfile(here, ['image.' exts{i}]);
    if exist(cand, 'file'), p = cand; return; end
end

% 2) any image in this folder
for i = 1:numel(exts)
    f = dir(fullfile(here, ['*.' exts{i}]));
    if ~isempty(f), p = fullfile(f(1).folder, f(1).name); return; end
end

% 3) synthesize a structure + texture test image
fprintf('No image found in folder - synthesizing a test image.\n');
[yy, xx] = meshgrid(linspace(0, 1, 384));
base = 0.55 + 0.40 * exp(-((xx - 0.62).^2 + (yy - 0.42).^2) / 0.008);
base = base + 0.30 * exp(-((xx - 0.30).^2) / 0.001) .* (yy > 0.2 & yy < 0.8);
tex  = 0.10 * sin(2*pi*xx*34) + 0.06 * sin(2*pi*yy*18) ...
     + 0.05 * sin(2*pi*(xx + yy)*48);
img  = repmat(min(1, max(0, base + tex)), 1, 1, 3);
p = fullfile(tempdir, 'image.png');
imwrite(img, p);
end
