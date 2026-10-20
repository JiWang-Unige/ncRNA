"""Launch genatator_distributed.py with one compatibility fix (run with cwd = pipeline repo).

Upstream bug (genatator-pipeline revision in MODEL_REVISIONS.json, transformers 4.57.6):
PretrainedConfig swallows the `dtype` attribute, so to_runtime_defaults() returns dtype=None
and ModelBank.load() rejects it. The pipeline only supports float32, so we restore
dtype='float32'. Verified 2026-10-20 that no other key of config.json changes on parsing.
"""
import os
import sys

sys.path.insert(0, os.getcwd())
import configuration_genatator_pipeline as cgp  # noqa: E402
import genatator_distributed  # noqa: E402

_orig = cgp.GenatatorPipelineConfig.to_runtime_defaults


def _patched(self):
    params = _orig(self)
    if params.get("dtype") in (None, "None"):
        params["dtype"] = "float32"
    return params


cgp.GenatatorPipelineConfig.to_runtime_defaults = _patched
genatator_distributed.GenatatorPipelineConfig.to_runtime_defaults = _patched

if __name__ == "__main__":
    sys.exit(genatator_distributed.main())
