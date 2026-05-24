# GLMCC container for reproducible hypothalamic FC analysis.
#
# Kobayashi et al. 2019 — Generalized Linear Model Cross-Correlation
# Repo: https://github.com/naokishibuya/GLMCC
#
# Build:
#   docker build -t glmcc:latest .
#
# Run (single animal):
#   docker run --rm -v $(pwd)/data:/data -v $(pwd)/results:/results \
#       glmcc --animal 171019 --condition lightON
#
# To pin a specific commit, uncomment the ARG line below and set:
#   docker build --build-arg GLMCC_COMMIT=<sha> -t glmcc:latest .

FROM python:3.12-slim

ARG GLMCC_COMMIT=main

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone https://github.com/naokishibuya/GLMCC.git && \
    cd GLMCC && \
    git checkout ${GLMCC_COMMIT}

# Install GLMCC Python dependencies (if requirements.txt exists)
WORKDIR /opt/GLMCC
RUN if [ -f requirements.txt ]; then pip install --no-cache-dir -r requirements.txt; fi

# Install additional deps for our pipeline
RUN pip install --no-cache-dir numpy scipy pandas

# Create data/results mount points
RUN mkdir -p /data/processed /results/glmcc

# Wrapper script
COPY scripts/run_glmcc_container.sh /usr/local/bin/run_glmcc
RUN chmod +x /usr/local/bin/run_glmcc

WORKDIR /workspace
ENTRYPOINT ["/usr/local/bin/run_glmcc"]
