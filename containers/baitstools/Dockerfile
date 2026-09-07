# BaitsTools -- https://github.com/campanam/BaitsTools
#
# No public conda package / biocontainer exists (checked 2026-09-07). BaitsTools
# ships only as a GitHub-hosted Ruby gem, so we build our own image.
#
#   docker build -t quay.io/sangerpathogens/baitstools:1.8.3 \
#       --build-arg BAITSTOOLS_VERSION=1.8.3 containers/baitstools/
#   docker push quay.io/sangerpathogens/baitstools:1.8.3
#
# gemspec requires Ruby >= 2.4.1; only runtime dep is `shell` 0.8.1 (GUI-only
# `tk` is skipped -- we install the non-gui executable set).

FROM ruby:3.2-slim

ARG BAITSTOOLS_VERSION=1.8.3

LABEL org.opencontainers.image.title="baitstools" \
      org.opencontainers.image.version="${BAITSTOOLS_VERSION}" \
      org.opencontainers.image.source="https://github.com/campanam/BaitsTools" \
      org.opencontainers.image.description="BaitsTools: hybridization capture bait design"

RUN apt-get update \
    && apt-get install -y --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 --branch "v${BAITSTOOLS_VERSION}" \
        https://github.com/campanam/BaitsTools /tmp/baitstools \
    && cd /tmp/baitstools \
    && gem build baitstools.gemspec \
    && gem install "./baitstools-${BAITSTOOLS_VERSION}.gem" \
    && cd / \
    && rm -rf /tmp/baitstools

# Smoke test: subcommand help must load.
RUN baitstools tilebaits --help > /dev/null

ENTRYPOINT []
CMD ["baitstools"]
