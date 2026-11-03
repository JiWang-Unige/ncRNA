#!/usr/bin/env python3
"""Run a TensorFlow script with its GPU memory capped:  tf_cap.py <MiB> <script.py> [args ...]
Sets a logical-device memory limit on GPU 0 before the script runs (Tiberius has no option for it)."""
import runpy
import sys

import tensorflow as tf

gpus = tf.config.list_physical_devices("GPU")
if gpus:
    tf.config.set_logical_device_configuration(gpus[0], [tf.config.LogicalDeviceConfiguration(memory_limit=int(sys.argv[1]))])
sys.argv = sys.argv[2:]
sys.path.insert(0, ".")
runpy.run_path(sys.argv[0], run_name="__main__")
