"""Download the GEO sample files listed in sources.csv (resumable; sizes checked).

    python data-raw/bayesprism/05_download.py <dest_dir> [--sizes]
"""
import csv, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

here = Path(__file__).resolve().parent
dest = Path(sys.argv[1])
only_sizes = "--sizes" in sys.argv


def url_of(gsm: str, fn: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{fn}"


def retry(f, tries=5):
    for i in range(tries):
        try:
            return f()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(10 * (i + 1))


def fetch(r):
    url = url_of(r["gsm"], r["file"])
    req = urllib.request.Request(url, method="HEAD")
    size = retry(lambda: int(urllib.request.urlopen(req, timeout=120).headers["Content-Length"]))
    out = dest / r["gse"] / r["file"]
    if only_sizes or (out.exists() and out.stat().st_size == size):
        return size
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    retry(lambda: urllib.request.urlretrieve(url, tmp))
    assert tmp.stat().st_size == size, url
    tmp.replace(out)
    print("got", r["model"], r["file"], size, flush=True)
    return size


rows = list(csv.DictReader(open(here / "sources.csv", encoding="utf-8")))
with ThreadPoolExecutor(6) as ex:
    total = sum(ex.map(fetch, rows))
print(f"total {total / 1e9:.2f} GB; all sizes match GEO")
