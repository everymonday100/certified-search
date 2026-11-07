#!/usr/bin/env python3
"""
Certified Search Verifier API (CS-025: API Hardening, CS-026: Rate Limiting).
HTTP-обёртка вокруг C++ верификатора с строгой валидацией, rate limiting и мониторингом.
"""
import json
import time
import hashlib
import subprocess
import tempfile
import os
from pathlib import Path
from typing import Optional, Literal
from datetime import datetime, timezone

from fastapi import FastAPI, HTTPException, Request, UploadFile, File, Form
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field, field_validator
import uvicorn

import sys
sys.path.append(str(Path(__file__).parent.parent))
from ledger.ledger_db import EvidenceLedger
from api.rate_limiter import SlidingWindowRateLimiter, MemoryMonitor

# --- Конфигурация ---
VERIFIER_BIN = "./core/verifier/build/krawczyk_lj13"
MAX_N_ATOMS = 55  # LJ55 — максимум по ТЗ
MAX_FILE_SIZE = 1024 * 1024  # 1 MB
MAX_COORD_VALUE = 1e6  # Защита от числового переполнения
ALLOWED_RADIUS_EXPONENTS = list(range(1, 13))  # 1e-1 ... 1e-12

LEDGER = EvidenceLedger(db_path="ledger/evidence.db")
RATE_LIMITER = SlidingWindowRateLimiter(max_requests=10, window_seconds=60)
MEMORY_LIMIT_MB = 1024.0  # 1 GB

app = FastAPI(
    title="Certified Search Verifier API",
    version="1.0.0",
    description="HTTP API для строгой верификации кластеров Леннарда-Джонса"
)

# --- Pydantic модели для валидации ---

class VerificationRequest(BaseModel):
    """JSON-запрос на верификацию."""
    coordinates: list[list[float]] = Field(..., description="Координаты атомов [[x,y,z], ...]")
    radius_exponent: int = Field(8, description="Радиус бокса = 10^(-radius_exponent)", ge=1, le=12)
    task_hash: Optional[str] = Field(None, description="Ожидаемый хэш задачи (для проверки соответствия)")
    
    @field_validator('coordinates')
    @classmethod
    def validate_coordinates(cls, v):
        n = len(v)
        if n < 2 or n > MAX_N_ATOMS:
            raise ValueError(f"Число атомов должно быть в диапазоне [2, {MAX_N_ATOMS}], получено {n}")
        for i, atom in enumerate(v):
            if len(atom) != 3:
                raise ValueError(f"Атом {i} должен иметь 3 координаты, получено {len(atom)}")
            for j, coord in enumerate(atom):
                if not (-MAX_COORD_VALUE <= coord <= MAX_COORD_VALUE):
                    raise ValueError(f"Атом {i}, координата {j}={coord} вне допустимого диапазона")
                if not (coord == coord):  # NaN check
                    raise ValueError(f"Атом {i}, координата {j} = NaN")
        return v
    
    @field_validator('radius_exponent')
    @classmethod
    def validate_radius(cls, v):
        if v not in ALLOWED_RADIUS_EXPONENTS:
            raise ValueError(f"radius_exponent должен быть в {ALLOWED_RADIUS_EXPONENTS}, получено {v}")
        return v

class VerificationResponse(BaseModel):
    status: Literal["CERTIFIED_STRICT_LOCAL_MINIMUM", "INCONCLUSIVE", "REJECTED"]
    energy_interval: Optional[list[float]] = None
    min_distance: Optional[float] = None
    kappa: Optional[float] = None
    pd_diagonal_lb: Optional[float] = None
    verification_time_sec: float
    ledger_record_hash: str

# --- Вспомогательные функции ---

def compute_file_hash(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()

def validate_xyz_content(content: str) -> tuple[int, list[list[float]]]:
    """Парсит XYZ-строку и валидирует её. Возвращает (N, coordinates)."""
    lines = content.strip().splitlines()
    if len(lines) < 2:
        raise ValueError("XYZ файл должен содержать минимум 2 строки (заголовок + комментарий)")
    
    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(f"Первая строка должна быть числом атомов, получено: {lines[0]!r}")
    
    if n_atoms < 2 or n_atoms > MAX_N_ATOMS:
        raise ValueError(f"Число атомов {n_atoms} вне диапазона [2, {MAX_N_ATOMS}]")
    
    if len(lines) < n_atoms + 2:
        raise ValueError(f"Ожидалось {n_atoms + 2} строк, получено {len(lines)}")
    
    coordinates = []
    for i in range(2, n_atoms + 2):
        parts = lines[i].split()
        if len(parts) < 4:
            raise ValueError(f"Строка атома {i-1} должна иметь формат 'Element x y z'")
        
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            raise ValueError(f"Строка атома {i-1} содержит нечисловые координаты: {parts[1:4]}")
        
        if not (x == x and y == y and z == z):
            raise ValueError(f"Атом {i-1} содержит NaN в координатах")
        if not (abs(x) < MAX_COORD_VALUE and abs(y) < MAX_COORD_VALUE and abs(z) < MAX_COORD_VALUE):
            raise ValueError(f"Атом {i-1} имеет координаты вне допустимого диапазона")
        
        coordinates.append([x, y, z])
    
    return n_atoms, coordinates

def run_verifier_binary(xyz_content: str, radius_exp: int) -> dict:
    """Запускает C++ верификатор и парсит вывод."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        temp_path = f.name
    
    try:
        start = time.perf_counter()
        result = subprocess.run(
            [VERIFIER_BIN, temp_path, str(radius_exp)],
            capture_output=True, text=True, timeout=60
        )
        elapsed = time.perf_counter() - start
        
        output = result.stdout
        error = result.stderr
        
        # Парсинг вывода
        status = "INCONCLUSIVE"
        energy_interval = None
        min_distance = None
        kappa = None
        pd_lb = None
        
        for line in output.splitlines():
            if "STATUS: CERTIFIED" in line:
                status = "CERTIFIED_STRICT_LOCAL_MINIMUM"
            elif "STATUS: REJECTED" in line:
                status = "REJECTED"
            elif "[energy] E in [" in line:
                # Парсим "[energy] E in [-44.3268, -44.3268]"
                try:
                    bounds = line.split("[")[2].split("]")[0].split(", ")
                    energy_interval = [float(bounds[0]), float(bounds[1])]
                except:
                    pass
            elif "[min_dist] lower bound=" in line:
                try:
                    min_distance = float(line.split("=")[1].strip())
                except:
                    pass
            elif "kappa_ub=" in line:
                try:
                    kappa_str = line.split("kappa_ub=")[1].split()[0]
                    kappa = float(kappa_str)
                except:
                    pass
            elif "min diagonal lower bound=" in line:
                try:
                    pd_lb = float(line.split("=")[1].strip())
                except:
                    pass
        
        return {
            "status": status,
            "energy_interval": energy_interval,
            "min_distance": min_distance,
            "kappa": kappa,
            "pd_diagonal_lb": pd_lb,
            "verification_time_sec": elapsed,
            "stdout": output,
            "stderr": error,
            "returncode": result.returncode
        }
    except subprocess.TimeoutExpired:
        return {
            "status": "REJECTED",
            "verification_time_sec": 60.0,
            "stdout": "",
            "stderr": "TIMEOUT",
            "returncode": -1
        }
    except Exception as e:
        return {
            "status": "REJECTED",
            "verification_time_sec": 0.0,
            "stdout": "",
            "stderr": str(e),
            "returncode": -1
        }
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)

# --- API Endpoints ---

@app.post("/verify", response_model=VerificationResponse)
async def verify_configuration(request: VerificationRequest):
    """
    Верификация конфигурации по JSON-запросу.
    """
    # Формируем XYZ-содержимое из координат
    n = len(request.coordinates)
    xyz_lines = [str(n), "API verification request"]
    for coord in request.coordinates:
        xyz_lines.append(f"X {coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f}")
    xyz_content = "\n".join(xyz_lines)
    
    config_hash = compute_file_hash(xyz_content.encode('utf-8'))
    task_hash = request.task_hash or f"sha256:api_request_{config_hash[:16]}"
    
    # Запуск верификатора
    result = run_verifier_binary(xyz_content, request.radius_exponent)
    
    # Запись в Ledger
    record_hash = LEDGER.add_attempt(
        task_hash=task_hash,
        config_hash=config_hash,
        status=result["status"],
        witness_data={
            "n_atoms": n,
            "radius_exponent": request.radius_exponent,
            "verification_time_sec": result["verification_time_sec"],
            "kappa": result.get("kappa"),
            "pd_diagonal_lb": result.get("pd_diagonal_lb")
        }
    )
    
    return VerificationResponse(
        status=result["status"],
        energy_interval=result.get("energy_interval"),
        min_distance=result.get("min_distance"),
        kappa=result.get("kappa"),
        pd_diagonal_lb=result.get("pd_diagonal_lb"),
        verification_time_sec=result["verification_time_sec"],
        ledger_record_hash=record_hash
    )

@app.post("/verify-xyz")
async def verify_xyz_file(file: UploadFile = File(...), radius_exponent: int = Form(8)):
    """
    Верификация XYZ-файла.
    """
    # Проверка размера файла
    content = await file.read()
    if len(content) > MAX_FILE_SIZE:
        raise HTTPException(status_code=413, detail=f"Файл слишком большой (максимум {MAX_FILE_SIZE} байт)")
    
    try:
        xyz_content = content.decode('utf-8')
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Файл должен быть в UTF-8")
    
    # Валидация XYZ
    try:
        n_atoms, coordinates = validate_xyz_content(xyz_content)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Невалидный XYZ: {str(e)}")
    
    # Проверка radius_exponent
    if radius_exponent not in ALLOWED_RADIUS_EXPONENTS:
        raise HTTPException(status_code=400, detail=f"radius_exponent должен быть в {ALLOWED_RADIUS_EXPONENTS}")
    
    config_hash = compute_file_hash(content)
    task_hash = f"sha256:api_xyz_{config_hash[:16]}"
    
    # Запуск верификатора
    result = run_verifier_binary(xyz_content, radius_exponent)
    
    # Запись в Ledger
    record_hash = LEDGER.add_attempt(
        task_hash=task_hash,
        config_hash=config_hash,
        status=result["status"],
        witness_data={
            "n_atoms": n_atoms,
            "radius_exponent": radius_exponent,
            "filename": file.filename,
            "verification_time_sec": result["verification_time_sec"]
        }
    )
    
    return {
        "status": result["status"],
        "energy_interval": result.get("energy_interval"),
        "min_distance": result.get("min_distance"),
        "kappa": result.get("kappa"),
        "pd_diagonal_lb": result.get("pd_diagonal_lb"),
        "verification_time_sec": result["verification_time_sec"],
        "ledger_record_hash": record_hash
    }

@app.get("/health")
async def health_check():
    """Проверка здоровья сервиса."""
    verifier_exists = Path(VERIFIER_BIN).exists()
    return {
        "status": "healthy" if verifier_exists else "degraded",
        "verifier_binary": VERIFIER_BIN,
        "verifier_exists": verifier_exists,
        "timestamp": datetime.now(timezone.utc).isoformat()
    }

@app.get("/ledger/verify-integrity")
async def verify_ledger_integrity():
    """Проверка целостности хэш-цепочки Ledger."""
    is_valid = LEDGER.verify_integrity()
    return {"integrity_valid": is_valid}

@app.get("/admin/rate-limit-stats")
async def rate_limit_stats():
    """Статистика rate limiter (для мониторинга)."""
    return {
        "stats": RATE_LIMITER.get_stats(),
        "memory_mb": round(MemoryMonitor.get_memory_mb(), 2),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

@app.post("/admin/rate-limit-reset")
async def rate_limit_reset(client_ip: Optional[str] = None):
    """Сброс rate limiter (для тестов)."""
    RATE_LIMITER.reset(client_ip)
    return {"status": "reset", "client_ip": client_ip}

# --- Middleware для rate limiting, memory monitoring и логирования ---

@app.middleware("http")
async def security_middleware(request: Request, call_next):
    start_time = time.perf_counter()
    
    # [1] Проверка памяти
    if not MemoryMonitor.check_memory_limit(MEMORY_LIMIT_MB):
        return JSONResponse(
            status_code=503,
            content={"detail": "Service temporarily unavailable: memory limit exceeded"},
            headers={"Retry-After": "60"}
        )
    
    # [2] Rate limiting (исключаем админские эндпоинты, чтобы можно было сбросить лимит)
    client_ip = request.client.host if request.client else "unknown"
    is_admin = request.url.path.startswith("/admin/")
    
    if not is_admin:
        allowed, info = RATE_LIMITER.is_allowed(client_ip)
        if not allowed:
            return JSONResponse(
                status_code=429,
                content={
                    "detail": "Rate limit exceeded",
                    "limit": info["limit"],
                    "retry_after": info["retry_after"],
                    "window_seconds": info["window_seconds"],
                },
                headers={
                    "Retry-After": str(info["retry_after"]),
                    "X-RateLimit-Limit": str(info["limit"]),
                    "X-RateLimit-Remaining": "0",
                }
            )
    else:
        # Для админских эндпоинтов создаем фиктивный info, чтобы код ниже не падал
        info = {"limit": 9999, "remaining": 9999}
    
    # Продолжаем обработку запроса
    response = await call_next(request)
    elapsed = time.perf_counter() - start_time
    
    # Добавляем rate limit заголовки к ответу
    response.headers["X-RateLimit-Limit"] = str(info["limit"])
    response.headers["X-RateLimit-Remaining"] = str(info["remaining"])
    
    # Структурированный лог
    log_entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "method": request.method,
        "path": request.url.path,
        "status_code": response.status_code,
        "elapsed_sec": round(elapsed, 4),
        "client_ip": client_ip,
        "rate_limit_remaining": info["remaining"],
        "memory_mb": round(MemoryMonitor.get_memory_mb(), 2),
    }
    print(json.dumps(log_entry), flush=True)
    
    return response

# --- Запуск ---

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000, log_level="info")