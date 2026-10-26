# OpenCode on top of the nightly's own image, for the CI steps where an AI model reads results.
#
# Same base as the nightly jobs: its layers are already on the TeamCity agents, and the model can
# use PHARE's Python stack (numpy, h5py, matplotlib, pyphare) as the tests do.
# Tag: registry.lpp.polytechnique.fr:5443/phare/opencode:<OPENCODE_VERSION>-<REVISION>. Bump REVISION
# when this file changes without a new OpenCode version: an existing tag is never rebuilt.

FROM registry.lpp.polytechnique.fr/phare/teamcity-fedora_dep:43

ARG OPENCODE_VERSION=1.18.32
ARG REVISION=1

# npm only to fetch the pinned OpenCode release (the same package as the GitHub workflow uses)
RUN dnf install -y --setopt=install_weak_deps=False nodejs npm && dnf clean all \
 && npm install --global "opencode-ai@${OPENCODE_VERSION}" && npm cache clean --force \
 && opencode --version

LABEL org.opencontainers.image.source="https://github.com/PHAREHUB/PHARE" \
      org.opencontainers.image.version="${OPENCODE_VERSION}-${REVISION}"
