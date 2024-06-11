"""
Generate random splits for partition settings
"""

from yaml import SafeDumper
import yaml
from jsonargparse import CLI
from typing import Union, Optional
import numpy as np

#   Output nones as empty strings


class NoAliasDumper(SafeDumper):
    def ignore_aliases(self, data):
        return True


NoAliasDumper.add_representer(
    type(None),
    lambda dumper, value: dumper.represent_scalar("tag:yaml.org,2002:null", ""),
)


def generate_split(
    system: str,
    T: Union[int, float],
    P: int,
    val_ratio: float = 0.2,
    test_ratio: float = 0.0,
    seed: int = 1234,
    stride: int = 1,
    batch_size: int = 64,
    subsample_random_seed: int = 42,
    max_n_frames: Optional[int] = None,
    max_epoch_samples: Optional[int] = None,
):
    """
    generate subsempling. If nax_n_frames is specified,
    use exprlicit indexing
    """

    h5_group = f"{system}_T={T:04d}K"  # h5 group name
    print(h5_group)

    molecules = [f"bead_{i}" for i in range(P)]  # list of replicas
    split_dict = {"val_ratio": val_ratio, "test_ratio": test_ratio, "seed": seed}
    if max_n_frames is None:
        detailed_indices = {molecule: split_dict.copy() for molecule in molecules}
    else:
        detailed_indices_val = dict()
        detailed_indices_train = dict()
        rng = np.random.default_rng(seed=seed)
        # Construct train metaset
        for molecule in molecules:
            frames = np.random.choice(max_n_frames, size=max_n_frames, replace=False)
            n_val_frames = int(val_ratio * max_n_frames)
            n_train_frames = max_n_frames - n_val_frames
            train_frames = list(np.copy(frames[:n_train_frames]))
            val_frames = list(np.copy(frames[n_train_frames:]))
            detailed_indices_val[molecule] = val_frames
            detailed_indices_train[molecule] = train_frames

    if max_n_frames is None:
        train = {
            "metasets": {
                h5_group: {
                    "molecules": molecules,
                    "detailed_indices": detailed_indices,
                    "stride": stride,
                }
            },
            "batch_sizes": {h5_group: batch_size},
            "subsample_random_seed": subsample_random_seed,
            "max_epoch_samples": max_epoch_samples,
        }

        full_partition = {"train": train, "val": train.copy()}
    else:
        train = {
            "metasets": {
                h5_group: {
                    "molecules": molecules,
                    "detailed_indices": detailed_indices_train,
                    "stride": stride,
                }
            },
            "batch_sizes": {h5_group: batch_size},
            "subsample_random_seed": subsample_random_seed,
            "max_epoch_samples": max_epoch_samples,
        }
        val = {
            "metasets": {
                h5_group: {
                    "molecules": molecules,
                    "detailed_indices": detailed_indices_train,
                    "stride": stride,
                }
            },
            "batch_sizes": {h5_group: batch_size},
            "subsample_random_seed": subsample_random_seed,
            "max_epoch_samples": max_epoch_samples,
        }

        full_partition = {"train": train, "val": val}
        print(full_partition)

    with open("partition_settings_global.yaml", "w") as f:
        yaml.dump(
            full_partition, f, Dumper=NoAliasDumper, indent=4, default_flow_style=False
        )
    return "partition_settings_global.yaml"


if __name__ == "__main__":
    CLI(generate_split)
