#!/usr/bin/env python3
"""M1 experiment 3 (docs/m1/PLAN_M1_AG.md): AlphaGenome fold-k models on held-out lncRNA / mRNA starts, in a Colab VM.

Started by ag_fold_run.sh, which uploads this file to /content/job.py, one checkpoint tar per model to
/content/<model>.tar and /content/job_args.json:
  {"mode": "smoke"|"full", "ag_commit": sha, "codes_sha256": {chrom: sha},
   "models": {"fold-k": {"tar_sha256": sha, "windows": {chrom: {k: [[row_id, offset, strand], ...]}}}}}
Install, sequence fetch (UCSC hg38, codes must equal the project's) and the fast prediction path follow
scripts/phase1/s2ag/ag_colab_job.py.  Per start, on the GPU, over [offset - 20 kb, offset + 20 kb] of the window and the
start's strand: for every CAGE and PRO-cap track the argmax (as offset from the start) and the maximum; the argmax of
the mean over CAGE tracks (primary score), of the max over CAGE tracks and of the mean over PRO-cap tracks.
Writes /content/out/<model>.npz + <model>.done (files with sha256) as each model finishes; tracks_<kind>_<strand>.csv.
Smoke: first model, first two windows; the fast path is checked against predict_sequence on one start.
"""
import hashlib
import json
import os
import subprocess
import sys
import time

T0 = time.time()
LOG = None
OUT = "/content/out"
W, STEP, H, HALF = 1 << 20, 1 << 19, 1 << 18, 20_000


def log(**kw):
    kw["t"] = round(time.time() - T0, 1)
    LOG.write(json.dumps(kw, default=str) + "\n")


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode:
        log(cmd=cmd, rc=r.returncode, stderr=r.stderr[-1500:])
        raise SystemExit(f"FAILED: {cmd}")
    return r.stdout


def gpu_ok():
    r = subprocess.run([sys.executable, "-c", "import jax; print(jax.devices()[0].platform)"], capture_output=True, text=True)
    return r.stdout.strip().endswith("gpu"), (r.stdout + r.stderr)[-600:]


def setup(args):
    url = f"git+https://github.com/google-deepmind/alphagenome_research.git@{args['ag_commit']}"
    pv = sh(f"{sys.executable} -c 'import google.protobuf as p; print(p.__version__)'").strip()
    sh(f"pip install -q '{url}'")
    sh(f"pip install -q hatchling 'protobuf=={pv}' grpcio-tools")
    sh(f"pip install -q --no-build-isolation --no-deps --force-reinstall '{url}'")
    sh(f"{sys.executable} -c 'import alphagenome_research.model.dna_model'")
    ok, msg = gpu_ok()
    if not ok:
        log(step="jax_no_gpu", msg=msg)
        sh("pip install -q -U 'jax[cuda12]'")
        ok, msg = gpu_ok()
    if not ok:
        raise SystemExit("FAILED: no GPU in jax")
    for m, spec in args["models"].items():
        got = sh(f"sha256sum /content/{m}.tar").split()[0]
        if got != spec["tar_sha256"]:
            raise SystemExit(f"FAILED: {m} checkpoint tar sha256 differs from ont")
        os.makedirs(f"/content/ckpt_{m}", exist_ok=True)
        sh(f"tar -xf /content/{m}.tar -C /content/ckpt_{m}")
    log(step="setup_done", nvidia=sh("nvidia-smi --query-gpu=name,memory.total --format=csv,noheader").strip(),
        versions=sh("pip list 2>/dev/null | grep -iE '^(jax|jaxlib|jax-cuda12-plugin|alphagenome|alphagenome.research|dm-haiku|orbax-checkpoint) '"))


_CODES = {}


def codes_of(chrom, want):
    import gzip
    import numpy as np
    if chrom in _CODES:
        return _CODES[chrom]
    path = f"/content/{chrom}.fa.gz"
    if not os.path.exists(path):
        sh(f"wget -q -O {path} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz")
    with gzip.open(path, "rb") as f:
        lines = f.read().split(b"\n")
    seq = b"".join(lines[1:])
    lut = np.full(256, 4, np.uint8)
    for i, ch in enumerate(b"ACGT"):
        lut[ch] = lut[ch + 32] = i
    codes = lut[np.frombuffer(seq, np.uint8)]
    got = hashlib.sha256(codes.tobytes()).hexdigest()
    log(step="codes", chrom=chrom, L=len(codes), match=got == want)
    if got != want:
        raise SystemExit(f"FAILED: {chrom} codes differ from the project's")
    _CODES.clear()
    _CODES[chrom] = codes
    return codes


class Runner:
    def __init__(self, model):
        import jax
        import jax.numpy as jnp
        import numpy as np
        from alphagenome.models.v1 import dna_model as api
        from alphagenome.models.v1 import dna_output
        from alphagenome_research.model import dna_model
        from alphagenome_research.model.metadata import metadata as metadata_lib
        self.jax, self.np, self.name = jax, np, model
        OT = dna_output.OutputType
        self.OT, self.org = OT, api.Organism.HOMO_SAPIENS
        self.req = (OT.CAGE, OT.PROCAP)
        t = time.time()
        self.model = dna_model.create(f"/content/ckpt_{model}", organism_settings={self.org: dna_model.OrganismSettings()})
        self.load_s = round(time.time() - t, 1)
        md = self.model._metadata[self.org]
        masks = metadata_lib.create_track_masks(md, requested_outputs=self.req, requested_ontologies=None)
        self.meta, idx = {}, {}
        for kind, ot in (("cage", OT.CAGE), ("procap", OT.PROCAP)):
            tab = md.get(ot)[masks[ot]].reset_index(drop=True)
            pad = np.flatnonzero(masks[ot])
            for s in "+-":
                sel = np.flatnonzero(tab["strand"].values == s)
                self.meta[(kind, s)] = tab.iloc[sel].reset_index(drop=True)
                idx[(kind, s)] = pad[sel]
        assert len(idx[("cage", "+")]) == len(idx[("cage", "-")]) and len(idx[("procap", "+")]) == len(idx[("procap", "-")]), \
            {k: len(v) for k, v in idx.items()}
        self.n_cage, self.n_pro = len(idx[("cage", "+")]), len(idx[("procap", "+")])
        dev = jax.local_devices()[0]
        self.dev = dev
        self.idx = {k: jax.device_put(v, dev) for k, v in idx.items()}
        self.oi = jax.device_put(np.full((1,), dna_model.convert_to_organism_index(self.org), np.int32), dev)
        self.neg = jax.device_put(np.asarray([False]), dev)
        self.sri = jax.device_put(md.strand_reindexing, dev)
        n = 2 * HALF + 1

        @jax.jit
        def stats(cage, pro, off, ci, pi):
            c = jax.lax.dynamic_slice_in_dim(cage, off - HALF, n, axis=0)[:, ci].astype(jnp.float32)
            p = jax.lax.dynamic_slice_in_dim(pro, off - HALF, n, axis=0)[:, pi].astype(jnp.float32)
            return (jnp.argmax(c, 0) - HALF, c.max(0), jnp.argmax(p, 0) - HALF, p.max(0),
                    jnp.argmax(c.mean(1)) - HALF, jnp.argmax(c.max(1)) - HALF, jnp.argmax(p.mean(1)) - HALF)
        self.stats = stats
        self.letters = np.frombuffer(b"ACGTN", np.uint8)

    def seq(self, codes, k):
        np = self.np
        s = k * STEP - H
        lo, hi = max(s, 0), min(s + W, len(codes))
        buf = np.full(W, 4, np.uint8)
        buf[lo - s:hi - s] = codes[lo:hi]
        return self.letters[buf].tobytes().decode()

    def predict(self, seq):
        m, np, jax = self.model, self.np, self.jax
        x = jax.device_put(np.asarray(m._one_hot_encoder.encode(seq))[np.newaxis], self.dev)
        p = m._predict(m._params, m._state, x, self.oi, requested_outputs=self.req,
                       negative_strand_mask=self.neg, strand_reindexing=self.sri)
        return p[self.OT.CAGE][0], p[self.OT.PROCAP][0]

    def starts(self, cage, pro, items):
        out = []
        for rid, off, s in items:
            r = self.stats(cage, pro, int(off), self.idx[("cage", s)], self.idx[("procap", s)])
            out.append((rid, [self.np.asarray(self.jax.device_get(x)) for x in r]))
        return out

    def public_check(self, seq, off, s):
        """argmax of the mean CAGE track (strand s) through predict_sequence vs the fast path, one start."""
        np = self.np
        o = self.model.predict_sequence(seq, organism=self.org, requested_outputs=self.req, ontology_terms=None)
        v = np.asarray(o.cage.values)[off - HALF:off + HALF + 1]
        st = o.cage.metadata["strand"].values
        pub = int(np.argmax(v[:, st == s].mean(1))) - HALF
        cage, pro = self.predict(seq)
        fast = self.starts(cage, pro, [["check", off, s]])[0][1]
        return {"public_mean_argmax": pub, "fast_mean_argmax": int(fast[4]), "n_public_tracks": int((st == s).sum()),
                "n_fast_tracks": self.n_cage}

    def peak_gb(self):
        try:
            st = self.dev.memory_stats() or {}
        except Exception:  # noqa: BLE001
            st = {}
        return round(st.get("peak_bytes_in_use", 0) / 2**30, 2)


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def done(tag, names, **extra):
    rec = {"files": {n: sha_file(f"{OUT}/{n}") for n in names}, **extra}
    json.dump(rec, open(f"{OUT}/{tag}.done.tmp", "w"), indent=1)
    os.replace(f"{OUT}/{tag}.done.tmp", f"{OUT}/{tag}.done")
    log(step="done", tag=tag, files=sorted(names))


def run_model(R, args, model, limit=None):
    np = R.np
    rid, cols = [], [[] for _ in range(7)]
    t, nwin = time.time(), 0
    for chrom, ws in args["models"][model]["windows"].items():
        if limit is not None and nwin >= limit:
            break
        codes = codes_of(chrom, args["codes_sha256"][chrom])
        for k, items in ws.items():
            if limit is not None and nwin >= limit:
                break
            cage, pro = R.predict(R.seq(codes, int(k)))
            for r, vals in R.starts(cage, pro, items):
                rid.append(r)
                for j, v in enumerate(vals):
                    cols[j].append(v)
            nwin += 1
            if nwin % 25 == 0:
                log(step="window", model=model, n=nwin, s=round(time.time() - t))
    names = ("cage_arg", "cage_max", "pro_arg", "pro_max", "agg_cage_mean_arg", "agg_cage_max_arg", "agg_pro_mean_arg")
    arr = {n: np.stack(c) for n, c in zip(names, cols)}
    for n in ("cage_max", "pro_max"):
        arr[n] = arr[n].astype(np.float16)
    for n in names:
        if n.endswith("arg"):
            arr[n] = arr[n].astype(np.int32)
    tag = model if limit is None else f"smoke_{model}"
    np.savez(f"{OUT}/{tag}.npz", rid=np.array(rid), **arr)
    return tag, nwin, round(time.time() - t, 1)


def main():
    global LOG
    LOG = open("/content/job.log", "a", buffering=1)
    args = json.load(open("/content/job_args.json"))
    log(step="start", mode=args["mode"], models=list(args["models"]))
    try:
        setup(args)
        os.makedirs(OUT, exist_ok=True)
        models = list(args["models"])[:1] if args["mode"] == "smoke" else list(args["models"])
        for model in models:
            if os.path.exists(f"{OUT}/{model}.done"):
                continue
            R = Runner(model)
            names = []
            if not os.path.exists(f"{OUT}/tracks_cage_+.csv"):
                for (kind, s), tab in R.meta.items():
                    tab.to_csv(f"{OUT}/tracks_{kind}_{s}.csv", index=False)
                names = [f"tracks_{k}_{s}.csv" for k in ("cage", "procap") for s in "+-"]
            log(step="model_loaded", model=model, load_s=R.load_s, n_cage=R.n_cage, n_procap=R.n_pro)
            if args["mode"] == "smoke":
                chrom, ws = next(iter(args["models"][model]["windows"].items()))
                codes = codes_of(chrom, args["codes_sha256"][chrom])
                k0, items = next(iter(ws.items()))
                t = time.time()
                chk = R.public_check(R.seq(codes, int(k0)), int(items[0][1]), items[0][2])
                tag, nwin, sec = run_model(R, args, model, limit=2)
                done("smoke", names + [f"{tag}.npz"], check=chk, public_s=round(time.time() - t - sec, 1), windows=nwin,
                     window_s=round(sec / max(nwin, 1), 2), peak_gb=R.peak_gb(), load_s=R.load_s,
                     n_cage=R.n_cage, n_procap=R.n_pro,
                     windows_total={m: sum(len(w) for w in v["windows"].values()) for m, v in args["models"].items()})
            else:
                tag, nwin, sec = run_model(R, args, model)
                done(model, names + [f"{tag}.npz"], windows=nwin, predict_s=sec, peak_gb=R.peak_gb(), load_s=R.load_s)
            del R
        log(step="end", total_s=round(time.time() - T0))
    except BaseException as e:  # noqa: BLE001 - the wrapper stops on this line
        log(step="FAILED", error=repr(e)[-1500:], total_s=round(time.time() - T0))
        raise


if __name__ == "__main__":
    main()
