"""Standalone interface to the research polarized compPSc solver.

All returned spectra use Q>0 parallel to the slab. Fortran COMMON blocks make
the native solver non-reentrant; calls are serialized within this module.
"""
import argparse
import ctypes as ct
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import tempfile
import numpy as np

SOURCE = Path(__file__).resolve().parent
SOURCE_FILES = ("compps.f", "mscomppsq.inc", "mscomppsj.inc", "standalone_logging.f")
LOCK = threading.RLock()

def build(directory="build", resolution=(51, 200, 7)):
    """Build a source-, compiler-, and mesh-specific shared library."""
    if len(resolution) != 3 or any(int(n) != n or n < 3 for n in resolution):
        raise ValueError("Resolution is (depth, energy, angle), each >=3")
    compiler = os.environ.get("COMPPSC_GFORTRAN") or shutil.which("gfortran")
    if not compiler:
        raise FileNotFoundError("Install gfortran or set COMPPSC_GFORTRAN")
    version = subprocess.check_output([compiler, "--version"], text=True)
    digest = hashlib.sha256((repr(resolution)+compiler+version+sys.platform).encode())
    for name in SOURCE_FILES:
        digest.update((SOURCE/name).read_bytes())
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    lib = directory/("libcomppsc_"+digest.hexdigest()[:16]+(".dylib" if sys.platform == "darwin" else ".so"))
    if lib.exists():
        return lib
    command = [compiler, "-O3", "-dynamiclib" if sys.platform == "darwin" else "-shared",
        "-fPIC", "-cpp", "-ffixed-line-length-none", "-std=legacy",
        "-fallow-argument-mismatch", "-fno-align-commons"]
    command += [f"-DCOMP_POL_MAX{k}={n}" for k,n in zip(("TAU","FRE","ANG"),resolution)]
    if sys.platform == "darwin":
        sdk = subprocess.check_output(["xcrun","--show-sdk-path"],text=True).strip()
        clang = Path(subprocess.check_output(["xcrun","--find","clang"],text=True).strip()).parent
        command += [f"-B{clang}/", "-isysroot",sdk,f"-Wl,-syslibroot,{sdk}"]
    fd, temporary = tempfile.mkstemp(prefix="comppsc-", suffix=lib.suffix, dir=directory)
    os.close(fd)
    command += ["compps.f","standalone_logging.f","-o",temporary]
    command += ["-framework","Accelerate"] if sys.platform == "darwin" else ["-lblas"]
    process = subprocess.run(command,cwd=SOURCE,capture_output=True,text=True)
    lib.with_suffix(".log").write_text(" ".join(command)+"\n"+process.stdout+process.stderr)
    if process.returncode:
        Path(temporary).unlink(missing_ok=True)
        process.check_returncode()
    os.replace(temporary, lib)
    return lib

def blackbody_bins(edges, kT_keV=.005):
    edges = np.asarray(edges, dtype=float)
    if (edges.ndim != 1 or len(edges)<2 or np.any(~np.isfinite(edges))
            or np.any(edges<=0) or np.any(np.diff(edges)<=0)):
        raise ValueError("Energy edges must be positive and strictly increasing")
    if not np.isfinite(kT_keV) or kT_keV <= 0:
        raise ValueError("Seed temperature must be positive")
    x,w = np.polynomial.legendre.leggauss(16)
    half = np.diff(edges)/2
    e = (edges[1:]+edges[:-1])[:,None]/2+half[:,None]*x
    ratio = e/kT_keV
    pdf = ratio**2/np.expm1(np.minimum(ratio,700))/(2*1.2020569031595943*kT_keV)
    pdf[ratio>=700] = 0
    return half*(pdf@w)

class PolarizedComppsc:
    def __init__(self, library, edges, seed_bin_probabilities, kTe_keV, tau):
        self.edges = np.asarray(edges,dtype=np.float32)
        self.seed = np.asarray(seed_bin_probabilities,dtype=np.float32)
        if (self.edges.ndim != 1 or not 3 <= len(self.edges) <= 5001 or np.any(~np.isfinite(self.edges))
            or np.any(self.edges<=0) or np.any(np.diff(self.edges)<=0)
            or self.seed.shape != (len(self.edges)-1,) or np.any(~np.isfinite(self.seed))
            or np.any(self.seed<0) or self.seed.sum()<=0 or not np.isfinite(kTe_keV+tau) or min(kTe_keV,tau)<=0):
            raise ValueError("Invalid physical parameters, edges, or seed probabilities")
        self.parameters = np.array([kTe_keV,2,-1,1000,tau,1,1,.5,1,0,1],dtype=np.float32)
        self.lib = ct.CDLL(str(Path(library).resolve()))
        self.ptr = ct.POINTER(ct.c_float)
        self.call = self.lib.xscomppspol2f_
        self.call.argtypes = [self.ptr,ct.POINTER(ct.c_int),self.ptr,ct.POINTER(ct.c_int)]+[self.ptr]*5
        self.call.restype = None
        self.stop = self.lib.comp_scan_stop_
        self.stop.argtypes = [ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
        self.stop.restype = None

    def spectrum(self, mu):
        """Return bin-integrated (I_up,Q_up,I_down,Q_down) per unit mu, audit."""
        if not 0 < mu <= 1:
            raise ValueError("mu must be in (0,1]")
        par = self.parameters.copy()
        par[7] = mu
        arrays = [self.seed.copy()]+[np.zeros_like(self.seed) for _ in range(4)]
        n,flag = ct.c_int(len(self.seed)),ct.c_int(0)
        residual,order,cap = ct.c_double(),ct.c_int(),ct.c_int()
        with LOCK:
            self.call(self.edges.ctypes.data_as(self.ptr),ct.byref(n),par.ctypes.data_as(self.ptr),
                ct.byref(flag),*(a.ctypes.data_as(self.ptr) for a in arrays))
            self.stop(ct.byref(residual),ct.byref(order),ct.byref(cap))
        audit = {"last_order":order.value,"order_cap":cap.value,
                 "last_relative_source_increment":residual.value}
        if not np.isfinite(residual.value) or residual.value>1e-7:
            raise ArithmeticError(f"Source iteration did not converge: {audit}")
        values = tuple(a.astype(float)*sign for a,sign in zip(arrays[:4],(1,-1,1,-1)))
        if any(np.any(~np.isfinite(a)) for a in values):
            raise ArithmeticError("Nonfinite spectrum")
        return values,audit

    def window(self, theta_min, theta_max, order=24):
        """Return differential photon spectra per keV per sr averaged over a window."""
        if not 0 <= theta_min < theta_max <= 90:
            raise ValueError("Require 0 <= theta_min < theta_max <= 90")
        lo,hi = np.cos(np.deg2rad([theta_max,theta_min]))
        nodes,weights = np.polynomial.legendre.leggauss(order)
        total = np.zeros((4,len(self.seed)))
        audits = []
        for x,w in zip(nodes,weights):
            values,audit = self.spectrum(lo+(x+1)*(hi-lo)/2)
            total += w*np.asarray(values)*(hi-lo)/2
            audits.append(audit)
        return total/(2*np.pi*(hi-lo)*np.diff(self.edges.astype(float))),audits

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kte",type=float,default=20.)
    p.add_argument("--tau",type=float,default=1.)
    p.add_argument("--seed-kt",type=float,default=.005)
    p.add_argument("--resolution",type=int,nargs=3,default=(51,200,7),metavar=("DEPTH","ENERGY","ANGLE"))
    p.add_argument("--mu",type=float,default=.5)
    p.add_argument("--output",type=Path,default=Path("output"))
    a = p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    edges = np.geomspace(1e-8,1e3,1761)
    solver = PolarizedComppsc(build(a.output/"build",tuple(a.resolution)),edges,
        blackbody_bins(edges,a.seed_kt),a.kte,a.tau)
    values,audit = solver.spectrum(a.mu)
    density = np.asarray(values)/(2*np.pi*np.diff(solver.edges.astype(float)))
    np.savetxt(a.output/"spectra.csv",np.column_stack((np.sqrt(solver.edges[:-1].astype(float)*solver.edges[1:]),density.T)),
        delimiter=",",header="E_keV,I_up,Q_up,I_down,Q_down",comments="")
    (a.output/"audit.json").write_text(json.dumps({"parameters":{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
        "iteration":audit,"seed_probability_on_input_grid":float(solver.seed.sum())},indent=2)+"\n")
    print(audit)

if __name__ == "__main__":
    main()
