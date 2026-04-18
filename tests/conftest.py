import sys, os
import glob
import re
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))


_REMAPPED_NAME_RE = re.compile(r"^(?P<prefix>.+)_remapped_(?P<assembly>[^/]+)\.csv$")


def pytest_addoption(parser):
    parser.addoption(
        "--results-dir",
        default=None,
        help="Root path of the pipeline output directory. "
             "Integration tests search the 'remapping/' and 'qc/' subdirectories automatically.",
    )
    parser.addoption(
        "--manifest",
        default=None,
        help="Path to the source Illumina manifest CSV that was fed to run_pipeline.sh. "
             "Required by integration tests; the pipeline does not persist this path in "
             "the results directory.",
    )


@pytest.fixture(scope="session")
def results_dir(request):
    d = request.config.getoption("--results-dir")
    if d is None:
        pytest.fail(
            "Integration tests require --results-dir. "
            "Run pytest with: --results-dir <path-to-results-folder>"
        )
    if not os.path.isdir(d):
        pytest.fail(f"--results-dir '{d}' does not exist or is not a directory.")
    return d


@pytest.fixture(scope="session")
def remapped_csv(results_dir):
    """Path to the single `{prefix}_remapped_{assembly}.csv` inside results_dir/remapping/."""
    pattern = os.path.join(results_dir, "remapping", "*_remapped_*.csv")
    matches = [m for m in glob.glob(pattern) if not m.endswith("_traced.csv")]
    if len(matches) == 0:
        pytest.fail(
            f"No '*_remapped_*.csv' found under {os.path.join(results_dir, 'remapping')}. "
            "Run run_pipeline.sh first."
        )
    if len(matches) > 1:
        pytest.fail(
            f"Ambiguous remapped CSV in {os.path.join(results_dir, 'remapping')}: "
            f"expected exactly one '*_remapped_*.csv', got {len(matches)}: {matches}"
        )
    return matches[0]


def _parse_remapped_name(path):
    m = _REMAPPED_NAME_RE.match(os.path.basename(path))
    if not m:
        pytest.fail(
            f"Remapped CSV filename {path!r} does not match "
            "expected '{prefix}_remapped_{assembly}.csv' pattern."
        )
    return m.group("prefix"), m.group("assembly")


@pytest.fixture(scope="session")
def assembly_label(remapped_csv):
    """Assembly tag embedded in the remapped CSV filename (e.g. 'equCab3', 'canFam3')."""
    _, assembly = _parse_remapped_name(remapped_csv)
    return assembly


@pytest.fixture(scope="session")
def manifest_path(request):
    """Absolute path to the source Illumina manifest CSV (from --manifest)."""
    p = request.config.getoption("--manifest")
    if p is None:
        pytest.fail(
            "Integration tests require --manifest. "
            "Run pytest with: --manifest <path-to-source-manifest.csv>"
        )
    if not os.path.isfile(p):
        pytest.fail(f"--manifest '{p}' does not exist or is not a file.")
    return os.path.abspath(p)
