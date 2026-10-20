#!/bin/bash
# Kaggle GPU benchmark of G1 log_partition at training scale (P1.0 performance check).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); d=$here/ncrna-g1-bench; mkdir -p $d
cat $here/../src/ncrna_annot/grammar/g1.py > $d/kernel.py
cat >> $d/kernel.py <<'PY'


if __name__ == "__main__":
    import json, time
    dev = "cuda"
    rng = np.random.default_rng(0)
    L = 24576
    seqs = ["".join(rng.choice(list("ACGT"), L)) for _ in range(8)]
    gates = torch.tensor(np.stack([gate_mask(encode(s)) for s in seqs]), device=dev)
    trans = G1Transitions().to(dev)
    out = {"gpu": torch.cuda.get_device_name(0), "torch": torch.__version__}
    for B in (1, 4, 8):
        for block in (64, 128, 256, 512):
            emit = torch.randn(B, L, S, device=dev, requires_grad=True)
            torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
            z = log_partition(emit, trans(), gates[:B], block=block); z.sum().backward(); torch.cuda.synchronize()
            t = time.time()
            for _ in range(3):
                emit.grad = None
                z = log_partition(emit, trans(), gates[:B], block=block); z.sum().backward()
            torch.cuda.synchronize()
            out[f"B{B}_block{block}"] = {"fwd_bwd_s": round((time.time() - t) / 3, 3),
                                         "peak_GB": round(torch.cuda.max_memory_allocated() / 1e9, 2),
                                         "finite": bool(torch.isfinite(z).all() and torch.isfinite(emit.grad).all())}
            print(f"B{B}_block{block}", out[f"B{B}_block{block}"], flush=True)
    json.dump(out, open("/kaggle/working/g1_bench.json", "w"), indent=1)
PY
cat > $d/kernel-metadata.json <<JS
{"id": "jiunige/ncrna-g1-bench", "title": "ncrna-g1-bench", "code_file": "kernel.py", "language": "python",
 "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false",
 "enable_internet": "false", "machine_shape": "", "dataset_sources": [], "competition_sources": [],
 "kernel_sources": [], "model_sources": []}
JS
echo $d
