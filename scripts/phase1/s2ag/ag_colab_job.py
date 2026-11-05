#!/usr/bin/env python3
"""S2-AG (PLAN_S2_AG.md): AlphaGenome all_folds on hg38 chromosomes, run inside a Colab GPU VM.

Started by ag_colab_run.sh, which uploads this file to /content/job.py, the checkpoint tar to /content/ag_ckpt.tar
and /content/job_args.json {"mode": "smoke"|"full", "chroms": [...], "codes_sha256": {chrom: sha}, "ag_commit": sha}.
Writes /content/job.log (progress lines; the last line carries "total_s") and /content/out/<tag>.done
({"files": {name: sha256}, ...}) once the files of <tag> are complete, so the wrapper can fetch them while the job runs.

Sequence: UCSC hg38 chromosome FASTA, encoded A0 C1 G2 T3 other 4; the code array must have the sha256 of the
project's codes (ont data/p1/*_codes.npy, data/splice_own/codes_chr7.npy), otherwise the job stops.
Windows: 1,048,576 bp, stride 524,288; each window keeps its central 524,288 bp, chromosome ends padded with N.
Kept per base (float16, no position offsets; ag_post.py picks them on ont):
  ag_<chrom>_splice.npy  [L, 4]  SPLICE_SITES classification, columns acc+, don+, acc-, don-
  ag_<chrom>_rna.npy     [L, 2]  RNA_SEQ, max over the + stranded tracks / over the - stranded tracks
The reduction runs on the GPU through the model's own jitted _predict (the public predict_sequence would move every
RNA-seq track to host memory); smoke mode checks it against predict_sequence on one window.
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
W = 1 << 20
STEP = W // 2
H = W // 4


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
    sh(f"pip install -q 'git+https://github.com/google-deepmind/alphagenome_research.git@{args['ag_commit']}'")
    ok, msg = gpu_ok()
    if not ok:
        log(step="jax_no_gpu", msg=msg)
        sh("pip install -q -U 'jax[cuda12]'")
        ok, msg = gpu_ok()
    if not ok:
        log(step="jax_no_gpu_after_fix", msg=msg)
        raise SystemExit("FAILED: no GPU in jax")
    os.makedirs("/content/ckpt", exist_ok=True)
    sh("tar -xf /content/ag_ckpt.tar -C /content/ckpt")
    log(step="setup_done", nvidia=sh("nvidia-smi --query-gpu=name,memory.total --format=csv,noheader").strip(),
        versions=sh("pip list 2>/dev/null | grep -iE '^(jax|jaxlib|jax-cuda12-plugin|alphagenome|alphagenome.research|dm-haiku|orbax-checkpoint) '"))


def codes_of(chrom, want):
    import numpy as np
    path = f"/content/{chrom}.fa.gz"
    if not os.path.exists(path):
        sh(f"wget -q -O {path} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz")
    import gzip
    with gzip.open(path, "rb") as f:
        lines = f.read().split(b"\n")
    assert lines[0].startswith(b">"), lines[0][:50]
    seq = b"".join(lines[1:])
    lut = np.full(256, 4, np.uint8)
    for i, ch in enumerate(b"ACGT"):
        lut[ch] = lut[ch + 32] = i
    codes = lut[np.frombuffer(seq, np.uint8)]
    got = hashlib.sha256(codes.tobytes()).hexdigest()
    log(step="codes", chrom=chrom, L=len(codes), sha256=got, match=got == want)
    if got != want:
        raise SystemExit(f"FAILED: {chrom} codes differ from the project's")
    return codes


class Runner:
    def __init__(self):
        import jax
        import jax.numpy as jnp
        import numpy as np
        from alphagenome.models.v1 import dna_model as api
        from alphagenome.models.v1 import dna_output
        from alphagenome_research.model import dna_model
        from alphagenome_research.model.metadata import metadata as metadata_lib
        self.jax, self.np = jax, np
        OT = dna_output.OutputType
        self.OT, self.org = OT, api.Organism.HOMO_SAPIENS
        self.req = (OT.SPLICE_SITES, OT.RNA_SEQ)
        t = time.time()
        self.model = dna_model.create("/content/ckpt", organism_settings={self.org: dna_model.OrganismSettings()})
        self.load_s = round(time.time() - t, 1)
        m = self.model
        self.md = m._metadata[self.org]
        masks = metadata_lib.create_track_masks(self.md, requested_outputs=self.req, requested_ontologies=None)
        smeta = self.md.get(OT.SPLICE_SITES)[masks[OT.SPLICE_SITES]].reset_index(drop=True)
        rmeta = self.md.get(OT.RNA_SEQ)[masks[OT.RNA_SEQ]].reset_index(drop=True)
        self.smeta, self.rmeta = smeta, rmeta
        log(step="tracks", splice=smeta.to_dict("records"), rna_cols=list(rmeta.columns), rna_n=len(rmeta),
            rna_strands=rmeta["strand"].value_counts().to_dict())
        spad = np.flatnonzero(masks[OT.SPLICE_SITES])
        rpad = np.flatnonzero(masks[OT.RNA_SEQ])
        cols = {}
        for i, r in smeta.iterrows():
            name = str(r["name"]).lower()
            kind = "acc" if "acceptor" in name else "don" if "donor" in name else None
            if kind and r["strand"] in "+-":
                cols[f"{kind}{r['strand']}"] = i
        order = ("acc+", "don+", "acc-", "don-")
        assert sorted(cols) == sorted(order), (cols, smeta.to_dict("records"))
        self.s_unpad = np.array([cols[c] for c in order])
        s_idx = spad[self.s_unpad]
        p_idx = rpad[np.flatnonzero(rmeta["strand"].values == "+")]
        m_idx = rpad[np.flatnonzero(rmeta["strand"].values == "-")]
        self.n_tracks = {"rna_plus": len(p_idx), "rna_minus": len(m_idx), "rna_unstranded": int((rmeta["strand"].values == ".").sum())}
        dev = jax.local_devices()[0]
        self.dev = dev
        self.consts = [jax.device_put(x, dev) for x in (s_idx, p_idx, m_idx)]
        self.oi = jax.device_put(np.full((1,), dna_model.convert_to_organism_index(self.org), np.int32), dev)
        self.neg = jax.device_put(np.asarray([False]), dev)
        self.sri = jax.device_put(self.md.strand_reindexing, dev)

        @jax.jit
        def reduce(sp, rn, si, pi, mi):
            s = sp[0][:, si].astype(jnp.float32)
            r = jnp.stack([rn[0][:, pi].max(-1), rn[0][:, mi].max(-1)], -1).astype(jnp.float32)
            return s, r
        self.reduce = reduce
        self.letters = np.frombuffer(b"ACGTN", np.uint8)

    def seq(self, codes, k):
        np = self.np
        s = k * STEP - H
        lo, hi = max(s, 0), min(s + W, len(codes))
        buf = np.full(W, 4, np.uint8)
        buf[lo - s:hi - s] = codes[lo:hi]
        return self.letters[buf].tobytes().decode()

    def fast(self, seq):
        m, np, jax = self.model, self.np, self.jax
        x = jax.device_put(np.asarray(m._one_hot_encoder.encode(seq))[np.newaxis], self.dev)
        p = m._predict(m._params, m._state, x, self.oi, requested_outputs=self.req,
                       negative_strand_mask=self.neg, strand_reindexing=self.sri)
        s, r = self.reduce(p[self.OT.SPLICE_SITES], p[self.OT.RNA_SEQ], *self.consts)
        return np.asarray(jax.device_get(s)), np.asarray(jax.device_get(r))

    def public(self, seq):
        np = self.np
        o = self.model.predict_sequence(seq, organism=self.org, requested_outputs=self.req, ontology_terms=None)
        s = np.asarray(o.splice_sites.values)[:, self.s_unpad]
        st = o.rna_seq.metadata["strand"].values
        rv = np.asarray(o.rna_seq.values)
        r = np.stack([rv[:, st == "+"].max(1), rv[:, st == "-"].max(1)], 1)
        return s, r

    def peak_gb(self):
        st = self.dev.memory_stats() or {}
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


def smoke(R, args):
    import zlib
    np = R.np
    chrom = args["chroms"][0]
    codes = codes_of(chrom, args["codes_sha256"][chrom])
    ks = list(range(60, 60 + args.get("smoke_windows", 3)))
    t = time.time()
    R.fast(R.seq(codes, ks[0]))
    first_s = round(time.time() - t, 1)
    times, outs = [], []
    for k in ks:
        t = time.time()
        outs.append(R.fast(R.seq(codes, k)))
        times.append(round(time.time() - t, 2))
    peak_fast = R.peak_gb()
    fs, fr = outs[0]
    diff = {"rna_max_value": float(fr.max()), "splice_max_value": float(fs.max())}
    t = time.time()
    try:
        ps, pr = R.public(R.seq(codes, ks[0]))
        diff.update(splice_max_abs=float(np.abs(ps - fs).max()), rna_max_abs=float(np.abs(pr - fr).max()))
    except Exception as e:  # noqa: BLE001 - e.g. host/device OOM of the public path; reported, not fatal
        diff["public_error"] = repr(e)[-500:]
    public_s = round(time.time() - t, 1)
    s16 = np.concatenate([o[0][H:H + STEP] for o in outs]).astype(np.float16)
    r16 = np.concatenate([o[1][H:H + STEP] for o in outs]).astype(np.float16)
    t = time.time()
    zr = {"splice": round(len(zlib.compress(s16.tobytes(), 1)) / s16.nbytes, 3),
          "rna": round(len(zlib.compress(r16.tobytes(), 1)) / r16.nbytes, 3), "zlib1_s": round(time.time() - t, 2)}
    os.makedirs(OUT, exist_ok=True)
    np.save(f"{OUT}/smoke_{chrom}_splice.npy", s16)
    np.save(f"{OUT}/smoke_{chrom}_rna.npy", r16)
    R.rmeta.to_csv(f"{OUT}/rna_seq_tracks.csv", index=False)
    R.smeta.to_csv(f"{OUT}/splice_sites_tracks.csv", index=False)
    n_win = sum(-(-L // STEP) for L in (64444167, 80373285, 138394717, 159345973))
    est = {"windows_4chrom": n_win, "est_predict_min": round(n_win * float(np.mean(times)) / 60, 1)}
    done("smoke", [f"smoke_{chrom}_splice.npy", f"smoke_{chrom}_rna.npy", "rna_seq_tracks.csv", "splice_sites_tracks.csv"],
         windows=ks, first_call_s=first_s, window_s=times, public_s=public_s, peak_gb_fast=peak_fast,
         peak_gb_after_public=R.peak_gb(), load_s=R.load_s, n_tracks=R.n_tracks, fast_vs_public=diff,
         zlib_ratio=zr, **est)
    log(step="smoke_summary", first_call_s=first_s, window_s=times, public_s=public_s, peak_gb=peak_fast,
        n_tracks=R.n_tracks, fast_vs_public=diff, zlib_ratio=zr, **est)


def full(R, args):
    np = R.np
    os.makedirs(OUT, exist_ok=True)
    for chrom in args["chroms"]:
        if os.path.exists(f"{OUT}/{chrom}.done"):
            continue
        codes = codes_of(chrom, args["codes_sha256"][chrom])
        L = len(codes)
        sp = np.zeros((L, 4), np.float16)
        rn = np.zeros((L, 2), np.float16)
        nk = -(-L // STEP)
        t = time.time()
        for k in range(nk):
            s, r = R.fast(R.seq(codes, k))
            lo, hi = k * STEP, min((k + 1) * STEP, L)
            sp[lo:hi] = s[H:H + hi - lo]
            rn[lo:hi] = np.minimum(r[H:H + hi - lo], 65000.0)
            if k % 25 == 0:
                log(step="window", chrom=chrom, k=k, of=nk, s=round(time.time() - t))
        pred_s = round(time.time() - t)
        names = [f"ag_{chrom}_splice.npy", f"ag_{chrom}_rna.npy"]
        np.save(f"{OUT}/{names[0]}", sp)
        np.save(f"{OUT}/{names[1]}", rn)
        done(chrom, names, L=L, windows=nk, predict_s=pred_s, peak_gb=R.peak_gb(),
             rna_clipped=int((rn >= 65000).sum()), splice_cols=["acc+", "don+", "acc-", "don-"], n_tracks=R.n_tracks)


def main():
    global LOG
    LOG = open("/content/job.log", "a", buffering=1)
    args = json.load(open("/content/job_args.json"))
    log(step="start", args={k: v for k, v in args.items() if k != "codes_sha256"})
    try:
        setup(args)
        R = Runner()
        log(step="model_loaded", load_s=R.load_s, n_tracks=R.n_tracks, splice_tracks=R.smeta.to_dict("records"))
        (smoke if args["mode"] == "smoke" else full)(R, args)
        log(step="end", total_s=round(time.time() - T0))
    except BaseException as e:  # noqa: BLE001 - the wrapper stops on this line
        log(step="FAILED", error=repr(e)[-1500:], total_s=round(time.time() - T0))
        raise


if __name__ == "__main__":
    main()
