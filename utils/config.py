from typing import Optional
from collections import Iterable


class RDFConfig:
    """
    Data class for RDF configuration
    """

    def __init__(
        self,
        r_range: Optional[Iterable[float]] = None,
        bin_width: Optional[float] = 0.005,
        n_bins: Optional[int] = None,
        periodic: Optional[bool] = True,
        opt: Optional[bool] = True,
    ):
        self.r_range = r_range
        self.bin_width = bin_width
        self.n_bins = n_bins
        self.periodic = periodic
        self.opt = opt


class BulkH2OConfig:
    def __init__(
        self,
        n_molecules: int,
        d: float,
        rdf_params_oo: RDFConfig,
        rdf_params_oh: RDFConfig,
        rdf_params_hh: RDFConfig,
        stride: int = 1,
    ):
        self.n_molecules = n_molecules
        self.d = d
        self.rdf_params_oo = rdf_params_oo
        self.rdf_params_oh = rdf_params_oh
        self.rdf_params_hh = rdf_params_hh
        self.stride = stride
