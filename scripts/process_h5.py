import h5py
import re
import numpy as np
from jsonargparse import CLI
import logging
from log_utils import log_execution_info

logger = logging.getLogger(__name__)


def preliminary_process(P: int, filename: str):

    f = h5py.File(filename, "a")
    # For each temperature, try to create an alias that does not have a space
    print(f.keys())

    # Proccessing stage
    for key in f.keys():
        temp_group = f[key]
        print(temp_group)
        for bead_key in temp_group.keys():
            if re.match("bead_", bead_key) is not None:
                bead = temp_group[bead_key]
                print(bead)
                atomic_forces = bead.get("for")
                atomic_forces = np.array(atomic_forces)
                spring_forces = bead.get("fsp")
                spring_forces = np.array(spring_forces)
                total_forces = (spring_forces + atomic_forces) / P
                try:
                    bead.create_dataset("total_forces", data=total_forces)
                except ValueError:
                    print("Dataset already exists")
    f.close()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger()
    log_execution_info(logger)

    CLI(preliminary_process)
    logger.info("Execution has ended successfully")
