# ==============================================================================
# STAGE 1: C++ Builder
# ==============================================================================
FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Установка зависимостей для сборки
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libflint-dev \
    libmpfi-dev \
    libmpfr-dev \
    libgmp-dev \
    libomp-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Копируем всю директорию verifier
COPY core/verifier/ ./core/verifier/

# Сборка
RUN mkdir -p core/verifier/build && cd core/verifier/build && \
    cmake .. && \
    make -j$(nproc)

# ==============================================================================
# STAGE 2: Python Runtime
# ==============================================================================
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Установка Python 3, pip и runtime-зависимостей C++
# ВАЖНО: libgomp1 необходим для бинарников, скомпилированных GCC с флагом -fopenmp
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    libflint18t64 \
    libmpfr6 \
    libgmp10 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Создаем виртуальное окружение для изоляции Python-пакетов
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Установка Python-зависимостей
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Копирование исходного кода Python
COPY api/ ./api/
COPY ledger/ ./ledger/
COPY scripts/ ./scripts/
COPY reference/ ./reference/

# Копирование скомпилированных C++ бинарников из Stage 1
COPY --from=builder /build/core/verifier/build/krawczyk_lj13 ./core/verifier/build/
COPY --from=builder /build/core/verifier/build/krawczyk_lj38 ./core/verifier/build/
COPY --from=builder /build/core/verifier/build/krawczyk_lj38_batch ./core/verifier/build/
COPY --from=builder /build/core/verifier/build/krawczyk_lj38_batch_opt ./core/verifier/build/

# Создание директории для Ledger
RUN mkdir -p /app/ledger && chmod 777 /app/ledger

EXPOSE 8000

# Запуск API
CMD ["uvicorn", "api.verifier_api:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "2"]