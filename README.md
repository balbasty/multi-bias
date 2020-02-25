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
and one mean image (prefixed by 'm');

### Command line

If more flexibility is required, `multibias` can be used as a function:
```{matlab}
>> [b,m] = multibias(x, opt);
```
The input `x` can be either:
- A 5D array: the fourth dimension should correspond to channels (and 
  therefore to bias fields), while the fifth channel should correspond to
  different contrasts. Note that in this model, bias fields are shared 
  across contrasts and differ across channels.
- A 2d cell array of filenames. Here the first dimension corresponds to 
  channels and the second to contrasts. Each file should contain a 3D
  volume. A 5D array is built from this combination of 3D files.
- Empty. In this case, files can be selected through the GUI.

The output are
- `b` is a 4D array of estimated bias field. The fourth dimension  
  corresponds to channels.
- `m` is a 5D array of combined 'mean' images. The fourth dimension 
  corresponds to contrasts, while the fourth dimension is of size 1.
- If no output argument is asked (`>> multibias(x, opt);`), an output 
  folder is asked for through the GUI and the output bias and mean images
  are writen on disk.

Several options are available, as fields of the `opt` structure:
| Name | Range | Default | Description |
|------|-------|---------|-------------|
| itermax | int > 0 | 16 | Maximum number of iterations |
| tol | float > 0 | 1E-3 | Gain threshold for early stopping |
| lambda | float > 0 | 1E3 | Regularisation factor per channel |
| sigma | float > 0 | NaN | Noise standard deviation per channel |
| vs | float > 0 | NaN/1 | Voxel size. <br>By default, it is inferred from the input headers. <br>If it cannot be inferred, the default value is 1. |
| verbose | int >= 0 | 1 | Verbosity level: 0=quiet / 1=print / 2=plot / 3=plot more |
| map | bool | false | Memory map input data (slower but saves RAM) |
| threads | int >= 0 <br> 'matlab' <br> 'automatic' | 'matlab' | Number of threads used by both Matlab and SPM: <br>'matlab' = Matlab's current settings <br>'automatic = Matlab's automatic setting |
| mask | bool | false | Mask of voxels to discard: <br>false = use all all voxels <br>true  = estimate and mask out background <br>array  = mask out using provided background mask |
| folder | char | '' | Output folder. Do not write output files if empty. |

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
