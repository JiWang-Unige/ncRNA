#!/home/ji/envs/kaggle/bin/python
"""Resumable, streaming download of one Kaggle kernel output file (the SDK buffers whole files in RAM, no resume).

usage: kaggle_output_resume.py <account> <owner/kernel> <file_name> <out_path>
Auth is identical to ~/bin/kaggle-as (token file path in env; the token and signed URLs are never printed).
"""
import os
import sys
import time
from pathlib import Path

import requests

account, kernel, fname, out = sys.argv[1:5]
profile = Path.home() / ".kaggle" / "accounts" / account
for n in ("KAGGLE_USERNAME", "KAGGLE_KEY", "KAGGLE_TOKEN", "KAGGLE_AUTH_METHOD"):
    os.environ.pop(n, None)
os.environ["KAGGLE_CONFIG_DIR"] = str(profile)
os.environ["KAGGLE_API_TOKEN"] = str(profile / "access_token")
from kaggle import api  # noqa: E402
from kaggle.api.kaggle_api_extended import ApiListKernelSessionOutputRequest  # noqa: E402

owner, slug = kernel.split("/")


def signed_url():
    token = None
    while True:
        with api.build_kaggle_client() as k:
            r = ApiListKernelSessionOutputRequest()
            r.user_name, r.kernel_slug = owner, slug
            if token:
                r.page_token = token
            resp = k.kernels.kernels_api_client.list_kernel_session_output(r)
        for it in resp.files or []:
            if it.file_name == fname:
                return it.url
        token = resp.next_page_token
        if not token:
            raise SystemExit(f"{fname} not in outputs of {kernel}")


total, fails = None, 0
while fails < 30:
    have = os.path.getsize(out) if os.path.exists(out) else 0
    if total is not None and have >= total:
        break
    try:
        h = {"Range": f"bytes={have}-"} if have else {}
        with requests.get(signed_url(), headers=h, stream=True, timeout=120) as r:
            if have and r.status_code != 206:
                raise SystemExit(f"server ignored Range (status {r.status_code})")
            r.raise_for_status()
            if r.status_code == 206:
                total = int(r.headers["Content-Range"].split("/")[-1])
            else:
                total = int(r.headers["Content-Length"])
            with open(out, "ab") as f:
                for chunk in r.iter_content(8 << 20):
                    f.write(chunk)
        print(f"segment ok, {os.path.getsize(out)}/{total}", flush=True)
    except (requests.RequestException, OSError, KeyError, ValueError) as e:
        fails += 1
        print(f"retry {fails}: {type(e).__name__} at {os.path.getsize(out) if os.path.exists(out) else 0}/{total}", flush=True)
        time.sleep(min(60, 5 * fails))
size = os.path.getsize(out)
print("DONE" if total and size == total else "INCOMPLETE", size, total, flush=True)
