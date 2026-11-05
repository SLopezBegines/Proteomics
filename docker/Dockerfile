# syntax=docker/dockerfile:1
# Proteomics DEP pipeline image: R 4.5.3 + Bioconductor 3.22 + DEP.
# Build from the repository root:  docker build -f docker/Dockerfile -t proteomics-dep:r4.5.3 .

# Base pinned by digest (rocker/r-ver:4.5.3, Ubuntu 24.04 "noble").
FROM rocker/r-ver:4.5.3@sha256:f5a72ffeea7ed592b2a5990f87101beef8ef3cae0e4ff7327088e80368f0aef9

# Dated CRAN snapshot used to build renv.lock (Posit Package Manager, precompiled binaries for noble).
ARG CRAN_SNAPSHOT=2026-06-30
ARG RENV_VERSION=1.3.1
ARG QUARTO_VERSION=1.11.9
ARG QUARTO_SHA256_AMD64=9d037baa3b38ffac7217edca1d2e53760b0fab4d5749220740a32056e03a1530
ARG QUARTO_SHA256_ARM64=8a96da89026c8cba0fc2607f62d99836eae8c5b65465b5d476af16f54c10d191

LABEL org.opencontainers.image.title="proteomics-dep" \
      org.opencontainers.image.description="DEP-based label-free proteomics pipeline (R 4.5.3, Bioconductor 3.22)" \
      org.opencontainers.image.source="https://github.com/SLopezBegines/Proteomics" \
      org.opencontainers.image.licenses="GPL-3.0"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    TZ=Etc/UTC

# System libraries needed at build time (compilation) and at run time (PPM binaries link against them).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev libxslt1-dev libicu-dev \
        libcairo2-dev libxt-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libfontconfig1-dev \
        libpng-dev libtiff-dev libjpeg-dev \
        libglpk-dev libgmp-dev libgit2-dev libnetcdf-dev libuv1-dev libsodium-dev libpcre2-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        graphviz libgraphviz-dev \
        curl git ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*

# Point runtime installs at the same dated snapshot as the lockfile.
RUN sed -i "s#^options(repos = .*#options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/noble/${CRAN_SNAPSHOT}'), download.file.method = 'libcurl')#" \
        "${R_HOME}/etc/Rprofile.site" \
    && grep -q "${CRAN_SNAPSHOT}" "${R_HOME}/etc/Rprofile.site"

# renv at the version recorded in the lockfile. The current CRAN release lives in
# src/contrib/; older ones move to Archive/, so try both.
RUN Rscript -e " \
      v <- '${RENV_VERSION}'; base <- 'https://cloud.r-project.org/src/contrib'; \
      urls <- c(sprintf('%s/renv_%s.tar.gz', base, v), sprintf('%s/Archive/renv/renv_%s.tar.gz', base, v)); \
      ok <- FALSE; \
      for (u in urls) { ok <- tryCatch({ install.packages(u, repos = NULL, type = 'source'); TRUE }, \
                                       warning = function(w) FALSE, error = function(e) FALSE); if (ok) break }; \
      stopifnot(ok, as.character(packageVersion('renv')) == v)"

# Restore the library. Only the lockfile is copied here, so this layer is cached
# until renv.lock changes (code edits never trigger a re-install).
COPY renv.lock /tmp/renv.lock
RUN Rscript -e " \
      options(Ncpus = max(1L, parallel::detectCores()), timeout = 900, \
              renv.config.repos.override = c(P3M = 'https://p3m.dev/cran/__linux__/noble/${CRAN_SNAPSHOT}', \
                                             CRAN = 'https://cloud.r-project.org')); \
      renv::restore(lockfile = '/tmp/renv.lock', library = .Library.site[1], prompt = FALSE, clean = FALSE)" \
    && rm -rf /tmp/renv.lock /tmp/Rtmp* /root/.cache/R

# Quarto CLI (bundles pandoc), checksum-verified. Placed after the renv layer so that
# changing it never invalidates the package restore.
ARG TARGETARCH
RUN set -eux; \
    case "${TARGETARCH:-amd64}" in \
      amd64) sha="${QUARTO_SHA256_AMD64}" ;; \
      arm64) sha="${QUARTO_SHA256_ARM64}" ;; \
      *) echo "unsupported arch ${TARGETARCH}"; exit 1 ;; \
    esac; \
    f="quarto-${QUARTO_VERSION}-linux-${TARGETARCH:-amd64}.tar.gz"; \
    curl -fsSL -o "/tmp/${f}" "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/${f}"; \
    echo "${sha}  /tmp/${f}" | sha256sum -c -; \
    mkdir -p /opt/quarto; tar -xzf "/tmp/${f}" -C /opt/quarto --strip-components=1; rm -f "/tmp/${f}"; \
    ln -s /opt/quarto/bin/quarto /usr/local/bin/quarto; \
    ln -s "$(find /opt/quarto -type f -name pandoc | head -n1)" /usr/local/bin/pandoc; \
    quarto --version; pandoc --version | head -n1

# Build-time verification: fail the build if the environment is not what we claim.
RUN Rscript -e " \
      stopifnot(R.version\$major == '4', startsWith(R.version\$minor, '5.3')); \
      stopifnot(as.character(BiocManager::version()) == '3.22'); \
      suppressPackageStartupMessages(library(DEP)); \
      stopifnot(as.character(packageVersion('DEP')) == '1.32.0'); \
      stopifnot(rmarkdown::pandoc_available()); \
      cat('OK: R', R.version.string, '| Bioconductor', as.character(BiocManager::version()), '| DEP', as.character(packageVersion('DEP')), '\n')"

# Non-root user. UID 1000 matches the usual host user, which keeps bind-mounted outputs writable.
RUN useradd --create-home --uid 1000 --shell /bin/bash proteomics
COPY --chown=proteomics:proteomics code/ /opt/proteomics/code/
WORKDIR /home/proteomics
USER proteomics
CMD ["R", "--quiet"]
