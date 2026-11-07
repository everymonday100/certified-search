#!/usr/bin/env python3
"""
Rate Limiter (CS-026).
Простой sliding-window rate limiter без внешних зависимостей.
Ограничивает количество запросов с одного IP за скользящее окно времени.
"""
import time
import threading
from collections import defaultdict
from typing import Optional

class SlidingWindowRateLimiter:
    """
    Sliding window rate limiter.
    Для каждого IP хранит timestamps последних запросов.
    Если за последние `window_seconds` было больше `max_requests` запросов — блокируем.
    """
    def __init__(self, max_requests: int = 10, window_seconds: int = 60):
        self.max_requests = int(max_requests)
        self.window_seconds = int(window_seconds)
        self._requests: dict[str, list[float]] = defaultdict(list)
        self._lock = threading.Lock()
    
    def is_allowed(self, client_ip: str) -> tuple[bool, dict]:
        """
        Проверяет, разрешён ли запрос от client_ip.
        Возвращает (allowed, info) где info содержит метаданные для заголовков.
        """
        now = time.time()
        window_start = now - self.window_seconds
        
        with self._lock:
            # Удаляем старые timestamps за пределами окна
            self._requests[client_ip] = [
                ts for ts in self._requests[client_ip] if ts > window_start
            ]
            
            current_count = len(self._requests[client_ip])
            remaining = max(0, self.max_requests - current_count)
            
            if current_count >= self.max_requests:
                # Блокируем: возвращаем время до разблокировки
                oldest = min(self._requests[client_ip])
                retry_after = int(oldest + self.window_seconds - now) + 1
                return False, {
                    "limit": self.max_requests,
                    "remaining": 0,
                    "retry_after": retry_after,
                    "window_seconds": self.window_seconds,
                }
            
            # Разрешаем: добавляем текущий timestamp
            self._requests[client_ip].append(now)
            return True, {
                "limit": self.max_requests,
                "remaining": remaining - 1,
                "retry_after": None,
                "window_seconds": self.window_seconds,
            }
    
    def get_stats(self) -> dict:
        """Возвращает статистику по всем IP (для мониторинга)."""
        now = time.time()
        window_start = now - self.window_seconds
        
        with self._lock:
            stats = {}
            for ip, timestamps in self._requests.items():
                active = [ts for ts in timestamps if ts > window_start]
                if active:
                    stats[ip] = {
                        "requests_in_window": len(active),
                        "limit": self.max_requests,
                        "window_seconds": self.window_seconds,
                    }
            return stats
    
    def reset(self, client_ip: Optional[str] = None):
        """Сбрасывает счётчики (для тестов или администрирования)."""
        with self._lock:
            if client_ip:
                self._requests.pop(client_ip, None)
            else:
                self._requests.clear()


class MemoryMonitor:
    """Простой монитор использования памяти процессом (для алертов)."""
    @staticmethod
    def get_memory_mb() -> float:
        """Возвращает RSS памяти процесса в МБ."""
        try:
            import resource
            # ru_maxrss в Linux — в КБ
            usage = resource.getrusage(resource.RUSAGE_SELF)
            return usage.ru_maxrss / 1024.0
        except Exception:
            return 0.0
    
    @staticmethod
    def check_memory_limit(limit_mb: float = 1024.0) -> bool:
        """Возвращает True, если память в пределах лимита."""
        return MemoryMonitor.get_memory_mb() < limit_mb