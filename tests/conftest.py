import sys, os
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))


def pytest_addoption(parser):
    parser.addoption(
        "--results-dir",
        default=None,
        help="Root path of the pipeline output directory. "
             "Integration tests search the 'remapping/' and 'qc/' subdirectories automatically.",
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
