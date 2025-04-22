apptainer exec --contain --cleanenv --pwd /FastCaloSim \
  -B $(pwd):/FastCaloSim \
  -B /tmp \
  -B $HOME/.gitconfig:$HOME/.gitconfig \
  docker://fcsproj/fastcalosim-alma9:main \
  bash -c "export PRE_COMMIT_HOME=/tmp/pre-commit-cache && pre-commit install && bash"
