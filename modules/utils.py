import os
import logging
from time import time
from subprocess import Popen, run, PIPE, STDOUT


def run_cmd(cmd: str, tool: str, log_stdout: bool = False) -> None:
    """Runs a shell command with pipefail set.

    Args:
        cmd (str): shell command to run.
        tool (str): tools being called.
        log_stdout (bool, optional): if the stdout should also be passed to the logs. Defaults to False.
    """
    with Popen(
        ["/bin/bash", "-c", "set -o pipefail; " + cmd],
        text=True,
        stdout=PIPE if log_stdout else None,
        stderr=STDOUT if log_stdout else PIPE,
    ) as p:
        for line in p.stdout if log_stdout else p.stderr:
            if (line := line.strip()) != "":
                logging.info(f"({tool}) {line}")


def validate_file(fpath: str) -> None:
    """Validates that the given path leads to legitimate file.

    Args:
        fpath (str): path to the file.
    """
    assert os.path.isfile(fpath), f"'{fpath}' is not a valid file..."


def return_nice_time(start: float, mins: bool = False) -> float:
    """Generated the delta time between now and a start unix time stamp.

    Args:
        start (float): unix time for the start.
        mins (bool, optional): if the time should be shown in minutes. Defaults to False.

    Returns:
        float: minutes or seconds of the time delta.
    """
    return round((time() - start) / (60 if mins else 1), 2)
