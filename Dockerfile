# QE Studio as a public web app (Hugging Face Spaces, Render, any Docker host).
# Build:  docker build -t qe-studio .
# Run:    docker run -p 7860:7860 -e QE_STUDIO_PASSWORD=choose-one qe-studio
FROM condaforge/miniforge3:latest

RUN conda create -y -n qe -c conda-forge qe numpy scipy matplotlib ase flask plotly gunicorn \
    && conda clean -afy

# Hugging Face Spaces runs containers as uid 1000
RUN if id -u 1000 >/dev/null 2>&1; then \
        old=$(id -nu 1000); usermod -l user -d /home/user -m "$old"; groupmod -n user "$(id -gn 1000)" || true; \
    else useradd -m -u 1000 user; fi
USER user

ENV HOME=/home/user \
    PATH=/opt/conda/envs/qe/bin:$PATH \
    QE_BIN=/opt/conda/envs/qe/bin \
    OMP_NUM_THREADS=1 \
    OMPI_MCA_smsc=^cma \
    OMPI_MCA_btl_vader_single_copy_mechanism=none \
    QE_STUDIO_PUBLIC=1 \
    QE_STUDIO_HOST=0.0.0.0 \
    QE_STUDIO_PORT=7860

WORKDIR /home/user/app
COPY --chown=user qe_studio/ .

EXPOSE 7860
# one worker: the job queue lives in this process; threads serve the UI while a job runs
CMD ["gunicorn", "-w", "1", "--threads", "8", "-b", "0.0.0.0:7860", "--timeout", "120", "app:app"]
