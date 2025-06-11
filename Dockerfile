FROM mambaorg/micromamba:1.5.6

# Add a dummy user entry to avoid `passwd` errors in Docker
RUN echo "micromamba:x:1000:1000::/home/micromamba:/bin/bash" >> /etc/passwd

COPY environment.yml /tmp/environment.yml
RUN micromamba create -n env -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/envs/env/bin:$PATH
