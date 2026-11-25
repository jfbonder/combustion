"""Helpers shared by the experiment scripts: provenance of the stored results
and process pools that carry the reaction rate to the workers.

Kept free of matplotlib so that every script (and the tests) can import it.
"""

import json
import pickle
import platform
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy

from combustion1d import model

NUMERICS = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def _git(*args):
    try:
        out = subprocess.run(("git",) + args, cwd=NUMERICS, capture_output=True, text=True,
                             timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def provenance():
    """Versions and commit that produced a result:
    {"git": <git rev-parse HEAD or None>, "git_dirty": <True if tracked files
    under numerics/ differ from that commit, None if unknown>, "python",
    "numpy", "scipy", "platform"}.

    No time stamp, so that deterministic results stay identical when they are
    regenerated with the same code and libraries.
    """
    status = _git("status", "--porcelain", "--untracked-files=no", "--", ".")
    return {
        "git": _git("rev-parse", "HEAD"),
        "git_dirty": None if status is None else bool(status),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": f"{platform.system()} {platform.machine()}",
    }


def with_provenance(obj, prov=None):
    """Copy of a result with the provenance attached, without changing its
    layout: a dict gets the key "provenance"; in a list of records (dicts)
    every record that does not have one gets it (records kept from an earlier
    run keep their own, or none if they predate this helper)."""
    prov = provenance() if prov is None else prov
    if isinstance(obj, dict):
        return obj if "provenance" in obj else dict(obj, provenance=prov)
    if isinstance(obj, (list, tuple)):
        return [dict(r, provenance=prov) if isinstance(r, dict) and "provenance" not in r else r
                for r in obj]
    return obj


def to_jsonable(o):
    """`default` hook of json.dumps for numpy values."""
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    raise TypeError(type(o))


def write_json(path, obj, indent=1):
    """Write obj (with its provenance, see with_provenance) as JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(with_provenance(obj), indent=indent, default=to_jsonable))


# ---------------------------------------------------------------------------
# Process pools
# ---------------------------------------------------------------------------

def worker_pool(n, reaction="inherit"):
    """ProcessPoolExecutor(n) whose workers use the same reaction rate as the
    parent process, whatever the start method.

    model.set_reaction changes module state of one process; workers started
    with "spawn" (macOS, Windows) re-import the package and would silently
    use the default f. Here the reaction (by default the one installed in the
    parent, model.get_reaction()) is installed in every worker by the pool
    initializer. A non-default reaction must be picklable (functions defined
    at module level); otherwise this raises a ValueError instead of running
    the workers with the wrong rate.
    """
    if reaction == "inherit":
        reaction = model.get_reaction()
    if reaction is not None:
        try:
            pickle.dumps(reaction)
        except Exception as ex:
            raise ValueError(
                "the installed reaction cannot be sent to worker processes (define f, df, M as "
                "module-level functions, not lambdas): " + repr(ex)) from ex
    return ProcessPoolExecutor(n, initializer=model.install_reaction, initargs=(reaction,))


if __name__ == "__main__":
    json.dump(provenance(), sys.stdout, indent=1)
    print()
