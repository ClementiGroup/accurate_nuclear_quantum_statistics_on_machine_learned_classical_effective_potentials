"""
File provides logging capabilities,
especially regarding the state of the git repository

Capabilities:
    Logging execution timestamp, location of the original script
    (including files pointed to with symbolic links),

"""

import logging
import inspect
import os
import sys
import subprocess


logger = logging.getLogger(__name__)


def log_execution_info(
    logger,
    log_python_info=True,
    log_repo_version=True,
):

    # Will craft several log messages and log them

    logger.info("** EXECUTION BEGINS ** ")

    # Add location of the script that calls the current function
    # (i.e. the script that is being executed)
    stack = inspect.stack()
    path_parent = os.path.abspath(stack[1].filename)
    source_code_loc = os.path.dirname(path_parent)
    message = f"Script location: {source_code_loc}"
    logger.info(message)
    # Log current working directory
    cwd = os.getcwd()
    message = f"Current working directory: {cwd}"
    logger.info(message)

    # Log the version of the repository (if it is a git repository)
    if log_repo_version:
        try:
            label = subprocess.run(
                [
                    "git",
                    "-C",
                    f"{source_code_loc}/.",
                    "describe",
                    "--always",
                    "--first-parent",
                    "--long",
                    "--abbrev=14",
                ],
                capture_output=True,
            ).stdout.decode(sys.stdout.encoding)
            logger.info(f"Current version of the git repository: {label}")
        except:
            logger.info("Was not able to get git repository info")

        # Log the status of the repository

        try:
            result = subprocess.run(
                ["git", "-C", source_code_loc, "status", "-s", "--porcelain"],
                capture_output=True,
            ).stdout
            if len(result) > 0:
                if result.split()[0] == "fatal:":
                    logger.info("No information about git tree state was found. \n")
                else:
                    logger.info("Repository tree state")
                    logger.info(result.decode(sys.stdout.encoding))
            else:
                logger.info("Clean working tree \n")
        except:
            logger.info("No information about git tree state was found. \n")


if __name__ == "__main__":
    # Test the function
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    log_execution_info(logger)
