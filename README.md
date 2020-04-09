# Multi-Bias

This code allows inferring multi-channel bias fields from, _e.g._, 
uncombined multi-coil MR images.

## Dependencies

This code is written in Matlab and depends on the 
[SPM12](https://www.fil.ion.ucl.ac.uk/spm/) software, which should be on
your Matlab path.

If you have access to the development version of SPM (_e.g._, if you work 
at the [FIL](https://www.fil.ion.ucl.ac.uk/spm/local/)), it is advised to 
recompile SPM with OpenMP activated. This can be done by following 
[these instructions](https://en.wikibooks.org/wiki/SPM), except that the 
option `USE_OPENMP=1` must be specified. For example, on linux:
```{shell}
cd /home/login/spm12/src
make distclean
make USE_OPENMP=1 && make install
make external-distclean
make external && make external-install
```

## Usage

### Graphical

Assuming that you have a series of uncombined coil images in nifti format, 
and that you like graphical interfaces, the most simple usage is:
```{matlab}
>> multibias;
```
You will first be asked to select all input files, then to choose an 
output directory. The algorithm will then run and write the results on 
disk. There will be as many bias field images as channels (prefixed by 'b')
and one mean image (prefixed by 'm').

### Command line

If more flexibility is required, `multibias` can be used as a function:
```{matlab}
>> [wb,m] = multibias(x, opt);
```
The input `x` can be either:
- A 5D numeric array: the fourth dimension should correspond to channels
  (and therefore to bias fields), while the fifth dimension should 
  correspond to different contrasts. Note that in this model, bias fields  
  are shared across contrasts and differ across channels.
- A 2d cell array of filenames (or file arrays or numeric arrays). 
  Here the first dimension corresponds to channels and the second to 
  contrasts. Each file should contain a 3D volume.
- Empty. In this case, files can be selected through the GUI.

The default output are
- `wb` is a cell of 3D arrays of estimated bias field, warped in native 
  (observed) space.
- `m` is a 5D array of combined 'mean' images. The fifth dimension 
  corresponds to contrasts, while the fourth dimension is of size 1.
- If no output argument is asked (`>> multibias(x, opt);`), an output 
  folder is asked for through the GUI and the output bias and mean images
  are writen on disk.

Several options are available, as fields of the `opt` structure:
| Name | Range | Default | Description |
|------|-------|---------|-------------|
| itermax | int > 0   | 16         | Maximum number of iterations |
| tol     | float > 0 | 1E-3       | Gain threshold for early stopping |
| armijo  | float > 0 | [8 4 2 1]  | Damping factors for Gauss-Newton |
| lambda  | float > 0 | 1E5        | Regularisation factor per channel |
| sigma   | float > 0 | NaN        | Noise standard deviation per channel |
| vs      | float > 0 | NaN/1      | Voxel size. <br>By default, it is inferred from the input headers. <br>If it cannot be inferred, the default value is 1. |
| verbose | int >= 0  | 1          | Verbosity level: 0=quiet / 1=print / 2=plot / 3=plot more |
| map     | bool      | false      | Memory map input data (slower but saves RAM) |
| threads | int >= 0 <br> 'matlab' <br> 'automatic' | 'matlab' | Number of threads used by both Matlab and SPM: <br>'matlab' = Matlab's current settings <br>'automatic = Matlab's automatic setting |
| mask    | bool      | false      | Mask of voxels to discard: <br>false = use all all voxels <br>true  = estimate and mask out background <br>array  = mask out using provided background mask |
| folder  | char      | ''         | Output folder. Do not write output files if empty. |
| output  | cellstr   | {'m' 'wb'} | List of returned arrays: <br> 'm':  mean image in mean space <br> 'wm': mean image warped to native spaces <br> 'b':  bias in mean space <br> 'wb': bias warped in native space <br> 'z':  log-bias in mean space <br> 'wz':  log-bias in native space <br> 'r':  world mean-to-native rigid mappings |
 
## Algorithm

This algorithm optimises a joint probability of the form `p(X | Y, Z) p(Z)`,
where:
- `X` contains `Nc x Ne` **observed** 3D volumes with contrast `e` and bias modulation `c`
- `Y` contains `Ne` **latent** demodulated 3D volumes
- `Z` contains `Nc` **latent** 3D log-bias fields
- `p(X{c,e} | Y{e}, Z{c}) = N(X{c,e} | exp(A{c,e}*Z{c}) .* A{c,e}*Y{e}, s{c,e}^2)`
- `A{c,e}` is a projection matrix that maps from a _mean_ space -- on which
  the latent volumes are defined -- to the observed space {c,e}
- `s{c,e}` is the noise standard eviation in the observe dimage {c,e}
- `p(Z{c}) = N(Z{c} | 0, inv(l{c}*L))`
- `L` is a precision matrix that captures the bending energy of the bias fields
- `l{c}` is a channel-specific regularisation factor

The mean space is obtained by computing the barycentre of the individual 
orientation matrices (after-coregistration).

The initial alignment is obtained by maximisation of the normalised 
mutual information between each observed image and the first one. 
To avoid biasing the orientation matrices towards the lattice of the first 
volume, all world-to-world alignement matrices are centred with respect 
to the mean of their Lie algebra representation:
> `R0   = exp(mean(Lie(R{:})))`
> `R{:} = R{:}/R0

The unmodulated images and log-bias fields are updated in turn using 
Gauss-Newton optimisation.

At each iteration, the log-bias fields are centered with-respect to their 
mean (weighted by regularisation factor).

## References

This algorithm has been described and used in a couple of ISMRM abstracts:

- **Estimation of net receive sensitivity - at 3T and 7T - for correction of inter-scan motion artefacts in R1 mapping.**  
[Yaël Balbastre](y.balbastre@ucl.ac.uk), [Nadège Corbin](n.corbin@ucl.ac.uk), [Martina F. Callaghan](m.callaghan@ucl.ac.uk)  
ISMRM 2020
- **Modelling RF coil-sensitivity induced artefacts in prospective motion-corrected accelerated 3D-EPI fMRI**  
[Nadine Graedel](n.graedel@ucl.ac.uk), [Yaël Balbastre](y.balbastre@ucl.ac.uk), [Nadège Corbin](n.corbin@ucl.ac.uk), [Oliver Josephs](o.josephs@ucl.ac.uk), [Martina F. Callaghan](m.callaghan@ucl.ac.uk)  
ISMRM 2020

## License

This software is released under the 
[GNU General Public License version 3](LICENSE) (GPL v3). As a result, 
you may copy, distribute and modify the software as long as you track 
changes/dates in source files. Any modifications to or software including 
(via compiler) GPL-licensed code must also be made available under the 
GPL along with build & install instructions.

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
