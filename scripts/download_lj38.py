#!/usr/bin/env python3
import urllib.request
import json
from pathlib import Path

def download_lj38():
    url = "https://www-wales.ch.cam.ac.uk/~wales/CCD/jon/structures/LJ/points/38?index=1"
    print(f"⬇️  Загрузка LJ38 глобального минимума из CCD...")
    
    req = urllib.request.Request(url, headers={"User-Agent": "certified-search/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = resp.read().decode('utf-8')
        
        numbers = [float(x) for x in data.split()]
        if len(numbers) != 114:  # 38 * 3
            raise ValueError(f"Ожидалось 114 координат, получено {len(numbers)}")
        
        out_dir = Path("reference/frozen/lj38")
        out_dir.mkdir(parents=True, exist_ok=True)
        
        xyz_file = out_dir / "lj38_global.xyz"
        with xyz_file.open("w") as f:
            f.write("38\n")
            f.write("LJ38 Global Minimum (Truncated Octahedron), E = -173.928427\n")
            for i in range(0, 114, 3):
                f.write(f"X {numbers[i]:.10f} {numbers[i+1]:.10f} {numbers[i+2]:.10f}\n")
        print(f"✅ Сохранено в {xyz_file}")
        
    except Exception as e:
        print(f"❌ Ошибка загрузки: {e}")

if __name__ == "__main__":
    download_lj38()