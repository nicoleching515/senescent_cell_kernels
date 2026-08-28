"""Resolve the DeepScence / senepy package DATA files without hardcoding one install.

Eight tracked scripts read

    /usr/local/lib/python3.11/dist-packages/DeepScence/data/coreGS_v2.csv
    /usr/local/lib/python3.11/dist-packages/senepy/data/

by absolute path.  That directory is the CONTAINER OVERLAY -- a second, undocumented
scientific stack, not the project venv, and precisely what master plan section 16.1 forbids
in bold; it is also what got wiped twice.  `pip install -r requirements.txt` into a venv does
not populate it at all, so on a clean rebuild those eight scripts fail (or, worse, read
whatever else is at that path).  See AUDIT_REPRODUCIBILITY D2.

Resolution order, first hit wins:
  1. an explicit override -- SASP_DEEPSCENCE_DATA / SASP_SENEPY_DATA;
  2. the package as importable by the RUNNING interpreter (importlib), i.e. the project venv
     when the drivers use code/_env.sh;
  3. the project venv's site-packages, even if this interpreter is not it;
  4. the container overlay, kept last so existing runs keep working unchanged.
Nothing is guessed silently: `where()` reports which of the four answered, and a miss raises
with all four candidates listed.

The two copies are currently byte-identical (coreGS_v2.csv md5 b981b9e9e730217d339306a709ada201
in both), so this changes no result today.  It changes whether a rebuilt clone can run at all.
"""
import importlib.util, os, sys

VENV = '/workspace/envs/sasp311/lib/python3.11/site-packages'
OVERLAY = '/usr/local/lib/python3.11/dist-packages'


def _candidates(pkg, env_var):
    out = []
    override = os.environ.get(env_var)
    if override:
        out.append(('$' + env_var, override))
    try:
        spec = importlib.util.find_spec(pkg)
    except (ImportError, ValueError):
        spec = None
    if spec is not None and spec.submodule_search_locations:
        out.append(('importlib (%s)' % sys.executable,
                    os.path.join(list(spec.submodule_search_locations)[0], 'data')))
    out.append(('project venv', os.path.join(VENV, pkg, 'data')))
    out.append(('container overlay', os.path.join(OVERLAY, pkg, 'data')))
    return out


def data_dir(pkg, env_var, must_contain=None):
    tried = []
    for how, d in _candidates(pkg, env_var):
        tried.append('%-34s %s' % (how, d))
        if os.path.isdir(d) and (must_contain is None
                                 or os.path.exists(os.path.join(d, must_contain))):
            return d, how
    raise SystemExit('cannot find the %s package data%s.  Tried:\n  %s\n'
                     'Install %s into the project environment '
                     '(/workspace/envs/sasp311) or set %s.'
                     % (pkg, '' if must_contain is None else ' (looking for %s)' % must_contain,
                        '\n  '.join(tried), pkg, env_var))


def deepscence_data(verbose=True):
    d, how = data_dir('DeepScence', 'SASP_DEEPSCENCE_DATA', 'coreGS_v2.csv')
    if verbose:
        print('DeepScence package data via %s: %s' % (how, d), file=sys.stderr)
    return d


def core_gs(verbose=True):
    """Path to DeepScence's coreGS_v2.csv."""
    return os.path.join(deepscence_data(verbose), 'coreGS_v2.csv')


def senepy_data(verbose=True):
    d, how = data_dir('senepy', 'SASP_SENEPY_DATA')
    if verbose:
        print('senepy package data via %s: %s' % (how, d), file=sys.stderr)
    return d + os.sep


def where():
    print('DeepScence coreGS_v2.csv ->', core_gs())
    print('senepy data/             ->', senepy_data())


if __name__ == '__main__':
    where()
