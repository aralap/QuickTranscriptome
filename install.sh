#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_DIR="${ROOT_DIR}/.miniconda"
ENV_NAME="quicktranscriptome"

if [[ ! -d "${CONDA_DIR}" ]]; then
  echo "Installing Miniconda into ${CONDA_DIR}..."
  UNAME="$(uname -s)"
  ARCH="$(uname -m)"

  if [[ "${UNAME}" == "Darwin" ]]; then
    if [[ "${ARCH}" == "arm64" ]]; then
      INSTALLER="Miniconda3-latest-MacOSX-arm64.sh"
    else
      INSTALLER="Miniconda3-latest-MacOSX-x86_64.sh"
    fi
  elif [[ "${UNAME}" == "Linux" ]]; then
    INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
  else
    echo "Unsupported OS: ${UNAME}" >&2
    exit 1
  fi

  curl -fsSL "https://repo.anaconda.com/miniconda/${INSTALLER}" -o "${ROOT_DIR}/miniconda.sh"
  bash "${ROOT_DIR}/miniconda.sh" -b -p "${CONDA_DIR}"
  rm -f "${ROOT_DIR}/miniconda.sh"
fi

source "${CONDA_DIR}/etc/profile.d/conda.sh"
if conda env list | awk '{print $1}' | awk -v env_name="${ENV_NAME}" '$1 == env_name {found=1} END {exit !found}'; then
  echo "Conda env ${ENV_NAME} already exists."
else
  echo "Creating Conda env ${ENV_NAME}..."
  conda env create -f "${ROOT_DIR}/environment.yml"
fi

echo
echo "Setup complete."
echo "Activate with:"
echo "  source \"${CONDA_DIR}/etc/profile.d/conda.sh\" && conda activate ${ENV_NAME}"
echo "Run pipeline with:"
echo "  python \"${ROOT_DIR}/quicktranscriptome.py\" run --reads-dir /path/to/fastqs"
