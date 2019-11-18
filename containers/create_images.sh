set -eux
sudo singularity build make_prg_dependencies.simg Singularity.make_prg_dependencies
sudo singularity build subsample.simg Singularity.subsample