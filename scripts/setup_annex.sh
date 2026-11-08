#!/usr/bin/env bash
# Attach OpenNeuro's annex metadata + public S3 remote to this fork so that
# `datalad get` / `git annex get` can fetch MRI content. Idempotent.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git remote get-url openneuro >/dev/null 2>&1 || \
  git remote add openneuro https://github.com/OpenNeuroDatasets/ds004725.git
git fetch openneuro git-annex:refs/remotes/openneuro/git-annex
git annex init "${HOSTNAME:-local}" >/dev/null
git config remote.origin.annex-ignore true
git annex enableremote s3-PUBLIC
git annex whereis sub-12000/func/sub-12000_task-rest_bold.nii
