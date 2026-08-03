"""HDF5 persistence helpers for FourierSolver.

This module is intentionally decoupled from the numerical solver: it defines
the :class:`HDF5Saver` mixin that adds an ``save_results_h5`` entry point
(plus a handful of small private helpers) to any class providing the right
attributes/methods. ``FourierSolver`` mixes it in via
``class FourierSolver(HDF5Saver):`` so the numerical core stays in
``fourier_solver.py`` and all I/O lives here.

Expected host attributes
------------------------
    n, L, dx, eps, dtype, model_params, grid_1d, gcalc

Expected host methods
---------------------
    solve_P_contribution(self) -> ndarray (N, N, N, N)
    solve_P_contribution__symmetry_optimized(self) -> ndarray (N, N, N, N)
    solve_Q_contribution(self) -> ndarray (N, N, N, N)
    solve_Q_contribution__symmetry_optimized(self) -> ndarray (N, N, N, N)

The host's ``gcalc`` is expected to expose ``calc_C_P_index`` and
``calc_C_Q_index`` with an ``out=`` keyword.
"""

import datetime as _dt
import os

import numpy as np

from .utils import log


class HDF5Saver:
    """Mixin: write requested fields to a self-describing HDF5 file.

    The single public entry point is :meth:`save_results_h5`. Each requested
    field is computed, written to disk, and then released before the next is
    computed, so peak RAM tracks the largest single field rather than their
    sum.

    Output paths can be either a concrete file or a directory; in the
    directory case :meth:`resolve_output_path` builds a deterministic,
    parameter-encoding filename via :meth:`default_filename`.
    """

    # Each save_* keyword maps to a short tag used in auto-generated
    # filenames. Iteration order also fixes the order of tags in the name.
    _SAVE_FLAG_TAGS = {
        "save_C_P":    "CP",
        "save_C_Q":    "CQ",
        "save_P":      "P",
        "save_P_diag": "Pd",
        "save_Q":      "Q",
        "save_Q_diag": "Qd",
    }

    # ------------------------------------------------------------------
    # Path / filename helpers
    # ------------------------------------------------------------------

    @classmethod
    def _save_tags_from_flags(cls, **save_flags):
        """Return the ordered list of short tags for True-valued ``save_*`` flags."""
        return [tag for flag, tag in cls._SAVE_FLAG_TAGS.items()
                if save_flags.get(flag, False)]

    def default_filename(self, *, method="symmetry", **save_flags):
        """Build a deterministic, parameter-encoding HDF5 filename.

        Floats use ``:g`` so trailing zeros disappear (1.0 -> "1") while
        irrationals like 7/3 keep enough precision to be unambiguous.
        Requires at least one True ``save_*`` flag so the filename
        unambiguously records what's inside.
        """
        save_tags = self._save_tags_from_flags(**save_flags)
        if not save_tags:
            raise ValueError(
                "default_filename: at least one save_* flag must be True."
            )
        parts = [
            "fourier",
            f"N{self.n}",
            f"L{self.L:g}",
        ]
        parts.extend(f"{name}={value:g}" for name, value in self.model_params.items())
        parts.extend([np.dtype(self.dtype).name, method])
        parts.append("save-" + "-".join(save_tags))
        return "_".join(parts) + ".h5"

    _H5_EXTENSIONS = (".h5", ".hdf5")

    def resolve_output_path(self, output, *, method="symmetry", **save_flags):
        """Resolve ``output`` to a concrete HDF5 file path.

        ``output`` is treated as a *file* iff its name ends in ``.h5`` or
        ``.hdf5``; otherwise it is treated as a *directory* and the
        auto-generated filename from :meth:`default_filename` is appended.
        This avoids ambiguity with :class:`pathlib.Path` (which strips
        trailing slashes) and with directories that don't yet exist.
        Accepts both ``str`` and :class:`os.PathLike`.
        """
        output = os.fspath(output)
        if output.lower().endswith(self._H5_EXTENSIONS):
            return output
        return os.path.join(
            output, self.default_filename(method=method, **save_flags)
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def save_results_h5(self, path, *,
                        save_C_P=False, save_C_Q=False,
                        save_P=False, save_P_diag=False,
                        save_Q=False, save_Q_diag=False,
                        method="symmetry"):
        """Compute the requested quantities and persist them to ``path``.

        Parameters
        ----------
        path : str or os.PathLike
            Destination HDF5 file. May also be a directory, in which case
            the filename is auto-generated via :meth:`resolve_output_path`.
        save_C_P, save_C_Q : bool
            Save all components of :math:`C_P` / :math:`C_Q` (no symmetry
            collapse) under the ``C_P/{ab}`` / ``C_Q/{abmn}`` groups.
        save_P, save_Q : bool
            Save the full 4-D pressure contribution as ``P_contrib`` /
            ``Q_contrib``.
        save_P_diag, save_Q_diag : bool
            Save the diagonal slice ``X[i, j, i, j]`` as a 2-D dataset.
        method : {"symmetry", "naive"}
            Which solver path to use for the P/Q contributions; ignored by
            the C_P/C_Q paths since those iterate every (a, b[, m, n])
            tuple anyway.

        Returns
        -------
        str
            Absolute path of the file that was written.

        Notes
        -----
        ``precompute()`` must have been called on the host first. The
        parent directory of ``path`` is created on demand.
        """
        # Local import keeps h5py optional for users who never call this.
        import h5py

        save_flags = {
            "save_C_P":    save_C_P,
            "save_C_Q":    save_C_Q,
            "save_P":      save_P,
            "save_P_diag": save_P_diag,
            "save_Q":      save_Q,
            "save_Q_diag": save_Q_diag,
        }
        if not any(save_flags.values()):
            raise ValueError(
                "save_results_h5: at least one save_* flag must be True."
            )
        if method not in {"symmetry", "naive"}:
            raise ValueError(
                f"save_results_h5: unknown method {method!r}; "
                f"expected 'symmetry' or 'naive'."
            )

        path = self.resolve_output_path(path, method=method, **save_flags)
        out_dir = os.path.dirname(os.path.abspath(path))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        with h5py.File(path, "w") as f:
            self._write_h5_attrs(f, method=method)
            f.create_dataset("grid_1d", data=self.grid_1d)

            # Cheap on-disk writes first: each component only needs a single
            # N^4 c_buffer in RAM (allocated inside the helper).
            if save_C_P:
                self._save_C_P_components(f)
            if save_C_Q:
                self._save_C_Q_components(f)

            # Heavy computations second: solve_*_contribution allocates a
            # complex Fourier-space buffer comparable in size to phi.
            if save_P or save_P_diag:
                P = self._compute_P_contribution(method)
                if save_P:
                    self._write_h5_dataset(f, "P_contrib", P)
                if save_P_diag:
                    self._write_h5_dataset(
                        f, "P_contrib_diag", self._diagonal_ijij(P)
                    )
                del P

            if save_Q or save_Q_diag:
                Q = self._compute_Q_contribution(method)
                if save_Q:
                    self._write_h5_dataset(f, "Q_contrib", Q)
                if save_Q_diag:
                    self._write_h5_dataset(
                        f, "Q_contrib_diag", self._diagonal_ijij(Q)
                    )
                del Q
        
        log(f"Saved {path}")
        return os.path.abspath(path)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _write_h5_attrs(self, f, *, method):
        f.attrs["N"]          = self.n
        f.attrs["L"]          = self.L
        f.attrs["dx"]         = self.dx
        f.attrs["eps"]        = self.eps
        f.attrs["dtype"]      = np.dtype(self.dtype).name
        f.attrs["method"]     = method
        f.attrs["created_at"] = _dt.datetime.now().isoformat(timespec="seconds")
        for name, value in self.model_params.items():
            f.attrs[name] = value

    @staticmethod
    def _write_h5_dataset(group, name, data):
        return group.create_dataset(name, data=data)

    def _save_C_P_components(self, f):
        """Save all 4 components of C_P as ``C_P/{a}{b}`` datasets.

        Each component is computed into a shared local buffer and then
        copied (synchronously, by h5py) into the file before the next
        component overwrites the buffer.
        """
        log("Saving C_P components", notime=True)
        g = f.create_group("C_P")
        c_buffer = np.empty((self.n,) * 4, dtype=self.dtype)
        for a in range(2):
            for b in range(2):
                log(f"Calculating C_P_index({a}, {b})")
                self.gcalc.calc_C_P_index(a, b, out=c_buffer)
                log(f"Writing C_P_index({a}, {b}) to h5")
                self._write_h5_dataset(g, f"{a}{b}", c_buffer)
                

    def _save_C_Q_components(self, f):
        """Save all 16 components of C_Q as ``C_Q/{a}{b}{m}{n}`` datasets."""
        log("Saving C_Q components", notime=True)
        g = f.create_group("C_Q")
        c_buffer = np.empty((self.n,) * 4, dtype=self.dtype)
        for a in range(2):
            for b in range(2):
                for m in range(2):
                    for n in range(2):
                        log(f"Calculating C_Q_index({a}, {b}, {m}, {n})")
                        self.gcalc.calc_C_Q_index(
                            a, b, m, n, out=c_buffer
                        )
                        self._write_h5_dataset(
                            g, f"{a}{b}{m}{n}", c_buffer
                        )

    def _compute_P_contribution(self, method):
        if method == "symmetry":
            return self.solve_P_contribution__symmetry_optimized()
        return self.solve_P_contribution()

    def _compute_Q_contribution(self, method):
        if method == "symmetry":
            return self.solve_Q_contribution__symmetry_optimized()
        return self.solve_Q_contribution()

    def _diagonal_ijij(self, x_rank4):
        """Return ``x[i, j, i, j]`` as a (N, N) array."""
        idx = np.arange(self.n)
        return x_rank4[idx[:, None], idx[None, :], idx[:, None], idx[None, :]]
