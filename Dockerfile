FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_PORT=7860 \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    SABR_MODEL_PATH=models/cas_subtype_extratrees.joblib

WORKDIR /app

RUN apt-get -o Acquire::Retries=5 update \
    && apt-get -o Acquire::Retries=5 install -y --no-install-recommends ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir --upgrade pip

COPY requirements.txt requirements-ml.txt requirements-external.txt ./
RUN python -m pip install --no-cache-dir \
    -r requirements.txt \
    -r requirements-ml.txt \
    -r requirements-external.txt

RUN useradd -m -u 1000 user

ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH

WORKDIR /home/user/app

COPY --chown=user app.py ./
COPY --chown=user crispr_phage_predictor/ ./crispr_phage_predictor/
COPY --chown=user assets/ ./assets/
COPY --chown=user data/examples/ ./data/examples/
COPY --chown=user scripts/ensure_model_artifact.py ./scripts/ensure_model_artifact.py

RUN mkdir -p models outputs/runs && chown -R user:user models outputs

USER user

EXPOSE 7860

CMD ["sh", "-c", "python scripts/ensure_model_artifact.py && streamlit run app.py --server.address \"$STREAMLIT_SERVER_ADDRESS\" --server.port \"$STREAMLIT_SERVER_PORT\""]
