.. _detailed-release-procedure:

Overview
========

The release procedure for ESMValTool is a fairly involved process (at the moment), so it
is important to be very well organized and to have documented each procedural steps, so that
the next release manager can follow said steps, and finalize the release without any delays.

The workflow below assumes an ESMValCore release candidate, or a completed stable release, have been released
and deployed on conda-forge and PyPI; it also assumes the release manager has access to accounts on LEVANTE [TODO insert info]

Procedural workflow
===================

1. Open an issue on GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~

First, open an issue on GitHub where the release workflow is documented. Name it something relevant like
"Recipe testing and comparison for release 2.x.x", and populate the isue description with information
about where the testing is taking place, what tools are used, and what versions, here is a template in Markdown:

# System parameters

Work done on DKRZ/Levante:

- submit files in `/home/b/b382109/submit`
- output in `/scratch/b/b382109/esmvaltool_output`

# System and Settings

## `conda`/`mamba`

```
(base) mamba --version
mamba 0.27.0
conda 22.9.0
```

## Git branch and state

Date: 25 October 2022 14:22 BST

```
(base) git status
On branch v2.7.x
Your branch is up to date with 'origin/v2.7.x'.
```

## Environment

On Levante:

```
mamba env create -n tool270Test -f environment.yml
conda activate tool270Test
```

## Environment file

Upload the full environment file here (use e.g. `conda env export > ToolEnv270Test.txt`).

[ToolEnv270Test.yml](https://github.com/ESMValGroup/ESMValTool/files/9860830/ToolEnv270Test.txt)

## Extraneous file movements

Mention any special data movement one needs to perfom to run the tests e.g.:

I moved the autoassess-specific files to `/home/b/b382109/autoassess_files` - run was succesful for AA recipes then :+1: 

## Modifications to configuration files

Example - I had to modify the `config-user.yml` file:

Added DKRZ downloaded data pool as:
```
  CMIP6:
    - /work/bd0854/DATA/ESMValTool2/CMIP6_DKRZ
    - /work/bd0854/DATA/ESMValTool2/download/CMIP6
  CMIP5:
    - /work/bd0854/DATA/ESMValTool2/CMIP5_DKRZ
    - /work/bd0854/b309141/additional_CMIP5
    - /work/bd0854/DATA/ESMValTool2/download/cmip5/output1
    - /work/bd0854/DATA/ESMValTool2/download/cmip5
```

as @schlunma and @remi-kazeroni have suggested :beer: 

# Recipe runs

Recipe runs results (as of **final** on 27 October 2022) are listed in https://github.com/ESMValGroup/ESMValTool/issues/2881#issuecomment-1291878142 (with very many thanks to @remi-kazeroni for running the impossible to run ones!) and are as follows:

- 122(121)*/127 successfully run recipes
- 0(1)*/127 failed with Diagnostic error, but fixed and rerun, but not yet PR-ed with the fix
- 2/127 that are missing data (for reals)
- 3/127 that have various issues (not missing data and not DiagnosticError)

`(*)` means not counting/counting the one that had a DiagnosticError but was fixed but not PR-ed

# Running the comparison

### Login and access to the DKRZ esmvaltool VM

Results from recipe runs are stored on the VM; login with:

```
ssh youraccount@esmvaltool.dkrz.de
```

### Get and install miniconda on VM

E.g. `scp Miniconda3-py39_4.12.0-Linux-x86_64.sh b382109@esmvaltool.dkrz.de:~` from a file already on Levante.

## Setting up the input files

If you wrote recipe runs output to Levante `/scratch` partition be aware that
the data will be removed after two weeks, so you will have to move the output data
to the `/work` partition, via e.g. a `nohup` job:

```
nohup cp -r /scratch/b/b382109/esmvaltool_output/* /work/bd0854/b382109/v270
```

`/work` is visible by the VM so you can run the compare tool straight on the VM.

**NOTE** do not store final release results on the VM including `/preproc/` dirs, the total
size for all the recipes output, including `/preproc/` dirs is in the 4.5TB ballpark,
much too high for the VM storage capacity

## Running compare tool at VM

- run date: 28 October 2022 (1st run)
- conda env: `tool270Compare`
- ESMValTool branch: `release270stable`
- prerquisite: `pip install imagehash`

### Input/output/run

- current: `/work/bd0854/b382109/v270` (contains `preproc/` dirs too, 122 recipes)
- reference: `/mnt/esmvaltool_disk2/shared/esmvaltool/v2.6.0rc4` (does not contain `preproc/` dirs)
- cmd: `nohup python ESMValTool/esmvaltool/utils/testing/regression/compare.py /mnt/esmvaltool_disk2/shared/esmvaltool/v2.6.0rc4 /work/bd0854/b382109/v270 > compare270output.txt`

Sanity check, as outputted by `compare.py`
```
Comparing recipe run(s) in:
/work/bd0854/b382109/v270
to reference in /mnt/esmvaltool_disk2/shared/esmvaltool/v2.6.0rc4
```
### First pass result

Running the `compare.py` results in a few recipes not-OK (NOK) wrt plots differing from previous release v2.6.0, summary in https://github.com/ESMValGroup/ESMValTool/issues/2881#issuecomment-1294735465

### Detailed plots inspection

Plots that differ for the 34 recipes that have them different is happening in https://github.com/ESMValGroup/ESMValTool/issues/2881#issuecomment-1295001054
