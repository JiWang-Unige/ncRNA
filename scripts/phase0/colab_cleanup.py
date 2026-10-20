"""Release Colab runtimes the CLI lost track of (e.g. `colab new` hit a ReadTimeout but the
server still assigned a VM, which then bills at the GPU rate).

  ~/.local/share/uv/tools/google-colab-cli/bin/python colab_cleanup.py          # list only
  ~/.local/share/uv/tools/google-colab-cli/bin/python colab_cleanup.py --orphans # unassign VMs not in local state
  ~/.local/share/uv/tools/google-colab-cli/bin/python colab_cleanup.py --all     # unassign every VM
"""
import sys

from colab_cli.common import state

client = state.client
known = {s.endpoint for s in state.store.list().values()}
assigned = client.list_assignments()
for a in assigned:
    orphan = a.endpoint not in known
    print(("ORPHAN " if orphan else "known  ") + a.endpoint)
    if "--all" in sys.argv or ("--orphans" in sys.argv and orphan):
        client.unassign(a.endpoint)
        print("  unassigned")
print("remaining", [a.endpoint for a in client.list_assignments()])
