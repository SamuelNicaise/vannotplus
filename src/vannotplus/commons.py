import logging as log
import os
from os.path import join as osj
import subprocess
import yaml

from vannotplus.family.ped9 import Ped


def set_log_level(verbosity):
    verbosity = verbosity.lower()
    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            f"Unknown verbosity level: {verbosity}\nPlease use any in: {configs.keys()}"
        )
    log.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


def load_config(config_file):
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)
    return config


def load_ped(config, app):
    ped_path = osj(config["ped_dir"], config["app_to_ped"][app])
    if not os.path.exists(ped_path):
        log.warning(f"No ped file found for app: {app}")
        ped = Ped()
    else:
        ped = Ped(ped_file=ped_path)
    return ped


def run_shell(cmd: str) -> None:
    """
    Run cmd in a separate shell
    Show stdout/stderr only if log level is log.DEBUG
    """
    log.debug(cmd)
    if log.root.level <= 10:
        redirect = None
    else:
        redirect = subprocess.DEVNULL
    subprocess.run(cmd, shell=True, stdout=redirect, stderr=redirect)
