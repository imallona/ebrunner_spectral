# Aim

Erich Brunner's and Giulia Moro's spectral project:
- probe designer
- probe designer app (under development)
- zUMIs analyzer

# Methods

## Probe designer

Lorem ipsum dolor sit amet

## Assemblies / GTF versions

Lorem ipsum dolor sit amet

## zUMIs

- Runs under `03_zumis/yaml/batch_1/no_barcodes_hm_0` are run without providing a barcodes' whitelist and with a Hamming Distance of 0 (i.e. no CB collapsing). In these runs, barcodes are assumed to be 9-spacer-9-spacer-9 followed by a 8nt-long UMI:

```
- BC(1-9,22-30,44-52)
- UMI(53-60)
```

# Repository structure

Developing at branch `develop` with occasional merges to `master`.

- `soft_installs.sh`: to set up the envs
- `01_proof_of_concept`: PoC for probe design (bash and Rmd)
- `02_browser`: genome browser trakcs generator / probe design (bash)
- `03_zumis`: utils to run zUMIs on the data, including the YAML files
    - `index_alien.sh`
    - `index_human.sh`
    - `index_alien.sh`
    - `index_human_mouse_alien_combined.sh`
    - `run_all.sh`
    - `yaml` folder including zUMIs yamls (refer to sherborne's paths)
        - `batch_1/no_barcodes_hm_0/human__20210923.B-o26015_1_3-SPECTRAL_unmod.yaml` etc
- `04_regex_harmonization`: utils to harmonize fastqs before running zUMIs

# How to run

Check `run_all.sh` and the software `soft_installs.sh` to get an idea of the setup (virtenvs, Rlibs etc).

```
./zUMIs.sh  -y test.yaml -d /home/ubuntu/sandbox/zUMIs/ 2>&1 | tee -a logs/output.log
```


# Notes/reminders

Currently being run in imlssherborne.uzh.ch.

Mind that that zUMIs and zUMIs dependencies are installed using `soft_installs.sh`.

# Started

Wed Aug 18 15:59:21 2021 +0200
