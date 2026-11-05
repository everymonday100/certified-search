FROM ubuntu:24.04
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential cmake \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY CMakeLists.txt ./
COPY src/ src/
COPY tests/ tests/
RUN cmake -B build -S . -DCMAKE_BUILD_TYPE=Release \
 && cmake --build build --parallel 2 \
 && ctest --test-dir build --output-on-failure
CMD ["bash"]
