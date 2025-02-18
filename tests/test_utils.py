import pathlib
import pytest
from time import time
from modules.utils import run_cmd, validate_file, return_nice_time


def test_run_cmd() -> None:
    test_cmd = "echo 'hello'"
    assert run_cmd(cmd=test_cmd, tool="", log_stdout=False) == None


def test_validate_file(tmp_path: pathlib.Path) -> None:
    fp = tmp_path / "test.txt"
    fp.touch()
    fp.write_text("hello text")
    assert validate_file(tmp_path / "test.txt") == None

    with pytest.raises(AssertionError, match=" is not a valid file..."):
        validate_file("fake_file.txt")


def test_nice_time() -> None:
    assert return_nice_time(start=time()) == pytest.approx(0.00)
    assert return_nice_time(start=time() - 10) == pytest.approx(10.00)
    assert return_nice_time(start=time() - 60, mins=True) == pytest.approx(1.00)
