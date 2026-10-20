# Log-Penalty-Based Multiscale Image Smoothing

Reference implementation of:

> Jianwu Long, Hao Luo, Shaoyi Wang. *Log-Penalty-Based Multiscale Image Smoothing for Edge Preservation and Texture Suppression.*

An edge-preserving smoothing method that removes textures and noise while keeping structural edges sharp.

## Files

| File | Description |
|---|---|
| `lp_smooth.m` | Core algorithm (self-contained single file) |
| `run_smoothing.m` | Minimal example: one image in, one smoothed image out |

## Requirements

- MATLAB R2018a or newer
- Image Processing Toolbox

## Quick Start

Put an image named **`image`** (e.g. `image.png` or `image.jpg`) in this folder, then run:

```matlab
run_smoothing
```

That's it — the smoothed result is shown next to the input and saved as `smoothed_image.png` in the same folder.

You can also smooth a specific file or use custom parameters:

```matlab
run_smoothing('my_photo.png')                    % any image path
run_smoothing([], 'lambda', 0.01, 'gamma', 3.0)  % custom parameters
S = lp_smooth('my_photo.png');                   % core function directly
imshow(S);
```

If no image is found in the folder, a synthetic structure+texture test image is generated automatically so the demo always runs.

Grayscale and RGB images are both supported. Input can be a file path or an image matrix (`uint8` or `double` in `[0,1]`); the output is `double` in `[0,1]`.

## Parameters

Only these four parameters are tunable; all other design choices are fixed to the configuration evaluated in the paper.

| Name | Default | Suggested range | Meaning |
|---|---|---|---|
| `lambda` | 0.01 | [0.005, 0.03] | Smoothing strength (most sensitive; larger = smoother) |
| `kappa` | 0.30 | [0.15, 0.40] | Gradient sensitivity of the edge-aware weight |
| `gamma` | 3.0 | [2.0, 5.0] | Texture gain of the wavelet prior |
| `rho` | 0.5 | [0.2, 0.7] | ADMM penalty parameter (keep ≤ 1) |
| `max_iter` | 20 | [10, 30] | ADMM iterations (converges within ~10) |

Example with custom parameters:

```matlab
run_smoothing('my_image.png', 'lambda', 0.02, 'gamma', 4.0);
```

## Citation

```bibtex
@article{long2026logpenalty,
  title   = {Log-Penalty-Based Multiscale Image Smoothing for Edge Preservation
             and Texture Suppression},
  author  = {Long, Jianwu and Luo, Hao and Wang, Shaoyi},
  journal = {(to appear)},
  year    = {2026}
}
```

## License

MIT License (see `LICENSE`).

## Contact

Hao Luo (corresponding author) — luohaocqut@stu.cqut.edu.cn
College of Computer Science and Engineering, Chongqing University of Technology
