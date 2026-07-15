#!/usr/bin/env python3
"""Python port of crumpling/lmp_swapstat4_main_bettercomments.cpp (core pipeline).

This is a faithful, top-to-bottom translation of the C++ ``main()`` covering the
*core* pipeline only:

  1. Read command line arguments
  2. Read the .ini parameter file
  3. Read the swap-file-names list
  4. Read the swap files (24 tab-separated columns per record)
  5. Detect subavalanches (cascades)
  6. Detect avalanches (mechanisms 0..5)
  7. Write the primary data outputs:
        <swapfile>__markup.txt
        <swapfile>__subavdata.txt
        <swapfile>__avdata.txt
        <outfname>__avdata.txt
        <outfname>_avstructure.txt

The histogram-based statistics sections of the C++ program (which depend on an
in-house math/stats library that is not part of this repository) are intentionally
NOT translated here.

Input format note: like the C++ program, the swap files must have exactly 24
tab-separated columns. The reduced ``*__swaps-short.txt`` files (6 columns) used by
the exploration notebook are NOT a valid input for this tool.

Usage (mirrors the original argv):
    python lmp_swapstat4_main.py <swapfnfile> <paramfile> <outfname>
"""

from __future__ import annotations

import math
import re
import sys
from typing import List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Small helpers mirroring v3_general / std behaviour
# ---------------------------------------------------------------------------

def _stod(s: str) -> float:
    """Mimic std::stod / strtod: parse the leading floating-point number."""
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        m = re.match(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?", s)
        return float(m.group(0)) if m else 0.0


def _stoi(s: str) -> int:
    """Mimic std::stoi / strtol: parse the leading integer."""
    s = s.strip()
    try:
        return int(s)
    except ValueError:
        m = re.match(r"[+-]?\d+", s)
        return int(m.group(0)) if m else 0


def fmt(value: float, prec: int) -> str:
    """Mimic ``std::fixed << std::setprecision(prec) << value`` for a double."""
    return f"{value:.{prec}f}"


def _max_pairwise_distance(pts: np.ndarray) -> float:
    """Max distance between any two of the given points (memory-safe, O(n^2))."""
    n = len(pts)
    if n < 2:
        return 0.0
    maxd2 = 0.0
    for i in range(n - 1):
        diff = pts[i + 1:] - pts[i]
        d2 = float(np.max(np.einsum("ij,ij->i", diff, diff)))
        if d2 > maxd2:
            maxd2 = d2
    return math.sqrt(maxd2)


# ---------------------------------------------------------------------------
# A tiny INI reader mimicking inih's INIReader (case-insensitive section+name)
# ---------------------------------------------------------------------------

class INIReader:
    def __init__(self, filename: str):
        self._values = {}
        self._error = 0
        try:
            with open(filename, "r") as f:
                lines = f.readlines()
        except OSError:
            self._error = -1
            return
        section = ""
        for raw in lines:
            line = raw.strip()
            if not line or line[0] in (";", "#"):
                continue
            if line[0] == "[":
                end = line.find("]")
                if end < 0:
                    self._error = 1
                    continue
                section = line[1:end].strip().lower()
                continue
            eq = line.find("=")
            if eq < 0:
                eq = line.find(":")
            if eq < 0:
                self._error = 1
                continue
            name = line[:eq].strip().lower()
            value = line[eq + 1:].strip()
            self._values[(section, name)] = value

    def ParseError(self) -> int:
        return self._error

    def _get(self, section: str, name: str):
        return self._values.get((section.lower(), name.lower()))

    def GetReal(self, section: str, name: str, default: float) -> float:
        v = self._get(section, name)
        return _stod(v) if v is not None else default

    def GetInteger(self, section: str, name: str, default: int) -> int:
        v = self._get(section, name)
        return _stoi(v) if v is not None else default

    def GetBoolean(self, section: str, name: str, default) -> bool:
        v = self._get(section, name)
        if v is None:
            return bool(default)
        return v.strip().lower() in ("true", "yes", "on", "1")


# ---------------------------------------------------------------------------
# Data structures (mirror the C++ structs)
# ---------------------------------------------------------------------------

class SwapRecord:
    __slots__ = (
        "nswap", "swapDirection", "atom1", "atom2",
        "bondCenter", "timestep", "time", "xsize", "ysize",
        "surfaceArea", "_likelihood_ijdtheta",
    )

    def __init__(self):
        self.nswap = 0
        self.swapDirection = 0
        self.atom1 = 0
        self.atom2 = 0
        self.bondCenter = np.zeros(3, dtype=float)
        self.timestep = 0.0
        self.time = 0.0
        self.xsize = 0.0
        self.ysize = 0.0
        self.surfaceArea = 0.0
        self._likelihood_ijdtheta = 0.0

    def set(self, ss: List[str]) -> int:
        if len(ss) != 24:
            return 1
        self.nswap = _stoi(ss[0])
        self.timestep = _stod(ss[2])
        self.time = _stod(ss[3])
        self.swapDirection = _stoi(ss[6])
        self.xsize = _stod(ss[7])
        self.ysize = _stod(ss[8])
        self.surfaceArea = _stod(ss[9])
        self.bondCenter = np.array([_stod(ss[17]), _stod(ss[18]), _stod(ss[19])], dtype=float)
        self.atom1 = _stoi(ss[22])
        self.atom2 = _stoi(ss[23])
        return 0

    def copy(self) -> "SwapRecord":
        r = SwapRecord()
        r.nswap = self.nswap
        r.swapDirection = self.swapDirection
        r.atom1 = self.atom1
        r.atom2 = self.atom2
        r.bondCenter = self.bondCenter.copy()
        r.timestep = self.timestep
        r.time = self.time
        r.xsize = self.xsize
        r.ysize = self.ysize
        r.surfaceArea = self.surfaceArea
        r._likelihood_ijdtheta = self._likelihood_ijdtheta
        return r


class Subavalanche:
    __slots__ = (
        "subavalancheID", "swaprecords_time_ID", "nswaps",
        "subavTimeP50", "subavTimeStart", "subavTimeEnd", "subavDuration",
        "subavCenterX", "subavCenterY", "subavRg", "subavI1", "subavI2",
        "subavDiameter",
    )

    def __init__(self):
        self.subavalancheID = 0
        self.swaprecords_time_ID: List[Tuple[float, int]] = []
        self.nswaps = 0
        self.subavTimeP50 = 0.0
        self.subavTimeStart = 0.0
        self.subavTimeEnd = 0.0
        self.subavDuration = 0.0
        self.subavCenterX = 0.0
        self.subavCenterY = 0.0
        self.subavRg = 0.0
        self.subavI1 = 0.0
        self.subavI2 = 0.0
        self.subavDiameter = 0.0


class Avalanche:
    __slots__ = (
        "avalancheID", "subavalanches_timeP50_ID", "subavalanches_timeStart_ID",
        "swaprecords_time_ID", "avNsubavs", "avNswaps", "avNswapsDir",
        "avTimeStart", "avTimeEnd", "avTimeAvg", "avThetaAvg", "avTimeP50",
        "avDuration", "avThetaDuration", "avTg", "avThetag", "avRg", "avDiameter",
    )

    def __init__(self):
        self.avalancheID = 0
        self.subavalanches_timeP50_ID: List[Tuple[float, int]] = []
        self.subavalanches_timeStart_ID: List[Tuple[float, int]] = []
        self.swaprecords_time_ID: List[Tuple[float, int]] = []
        self.avNsubavs = 0
        self.avNswaps = 0
        self.avNswapsDir = 0
        self.avTimeStart = 0.0
        self.avTimeEnd = 0.0
        self.avTimeAvg = 0.0
        self.avThetaAvg = 0.0
        self.avTimeP50 = 0.0
        self.avDuration = 0.0
        self.avThetaDuration = 0.0
        self.avTg = 0.0
        self.avThetag = 0.0
        self.avRg = 0.0
        self.avDiameter = 0.0


# ---------------------------------------------------------------------------
# Shared avalanche property computation (mirrors the duplicated C++ blocks).
#
# The C++ code recomputes these properties identically in mechanism 0 (inline,
# on a time-sorted swap list) and in mechanisms 1..5 (post-loop, on the
# append-ordered swap list). To preserve that subtlety, mechanism 0 sorts its
# swaprecords_time_ID before calling this; mechanisms 1..5 do not.
# ---------------------------------------------------------------------------

def compute_av_properties(av: Avalanche, swaprecords: List[SwapRecord], positions: np.ndarray):
    times = [t for (t, _) in av.swaprecords_time_ID]
    ids = [i for (_, i) in av.swaprecords_time_ID]
    n = len(times)

    av.avNsubavs = len(av.subavalanches_timeP50_ID)
    av.avNswaps = n
    av.avNswapsDir = -sum(swaprecords[i].swapDirection for i in ids)

    av.avTimeStart = times[0]
    av.avTimeEnd = times[-1]
    av.avDuration = av.avTimeEnd - av.avTimeStart
    av.avThetaDuration = math.log10(av.avTimeEnd) - math.log10(av.avTimeStart)

    logtimes = [math.log10(t) for t in times]
    av.avTimeAvg = sum(times) / n
    av.avThetaAvg = sum(logtimes) / n

    pts = positions[ids]
    avcenter = pts.mean(axis=0)

    if n % 2 == 1:
        av.avTimeP50 = times[n // 2]
    else:
        av.avTimeP50 = (times[n // 2] + times[n // 2 - 1]) * 0.5

    tg = 0.0
    thetag = 0.0
    for t in times:
        dt = av.avTimeAvg - t
        tg += dt * dt
    for lt in logtimes:
        dtheta = av.avThetaAvg - lt
        thetag += dtheta * dtheta
    dr = pts - avcenter
    rg = float(np.einsum("ij,ij->", dr, dr))

    av.avDiameter = _max_pairwise_distance(pts)
    av.avRg = math.sqrt(rg / n)
    av.avTg = math.sqrt(tg / n)
    av.avThetag = math.sqrt(thetag / n)


def main(argv: List[str]) -> int:
    # =====================================================================
    # Technical block for initialization and data parsing
    # =====================================================================

    # $ Read command line arguments
    if len(argv) < 3:
        print("\nNot enough arguments!", end="")
        return 1
    swapfnfile = argv[0]
    paramfname = argv[1]
    outfname = argv[2]
    # / Read command line arguments

    # $ Read ini file
    print("\nRead ini file...", end="")
    ini = INIReader(paramfname)
    if ini.ParseError() != 0:
        print(f'Error parsing parameters file "{paramfname}"')
        return 1

    # System
    system_posLeft_x = ini.GetReal("System", "system_posLeft_x", 0.00)
    system_posLeft_y = ini.GetReal("System", "system_posLeft_y", 0.00)

    # SubavalancheDetection
    avgBondLength = ini.GetReal("SubavalancheDetection", "avgBondLength", 0.00)
    subavalanche_dt_threshold = ini.GetReal("SubavalancheDetection", "subavalanche_dt_threshold", 10.00)
    subavalanche_dr_threshold = ini.GetReal("SubavalancheDetection", "subavalanche_dr_threshold", 100.00)

    # AvalancheDetection
    #  0 - by dt/t cutoff
    #  1 - by dt/t and r cutoffs, closest in space wins
    #  2 - by dt/t and r cutoffs, bigger wins
    #  3 - by dt/t and r cutoffs plus minimum dt of the avalanche
    #  4 - by dt/t and r/dt > v cutoffs, bigger wins
    #  5 - by dt/t and r/dt > v cutoffs, biggest choice
    avalancheDetectionMechanism = ini.GetInteger("AvalancheDetection", "avalancheDetectionMechanism", 0)
    avDetection_subavNswaps_threshold = ini.GetInteger("AvalancheDetection", "avDetection_subavNswaps_threshold", 1)
    avDetection_dtt_threshold = ini.GetReal("AvalancheDetection", "avDetection_dtt_threshold", 0.01)
    avDetection_dt_minimum = ini.GetReal("AvalancheDetection", "avDetection_dt_minimum", 0.01)
    avDetection_dr_threshold = ini.GetReal("AvalancheDetection", "avDetection_dr_threshold", 100.00)
    avDetection_velocity_threshold = ini.GetReal("AvalancheDetection", "avDetection_velocity_threshold", 20.00)

    # Statistics
    analysis_time_min = ini.GetReal("Statistics", "analysis_time_min", 0.00)
    analysis_time_max = ini.GetReal("Statistics", "analysis_time_max", -1.00)

    # AvalancheStat (only the avStructure_* parameters are used by the core pipeline)
    avStructure_nswaps_threshold = ini.GetInteger("AvalancheStat", "avStructure_nswaps_threshold", 0)
    avStructure_dt = ini.GetReal("AvalancheStat", "avStructure_dt", 1.0)
    avStructure_dlogt = ini.GetReal("AvalancheStat", "avStructure_dlogt", 1.0)
    avStructure_dlogtnorm = ini.GetReal("AvalancheStat", "avStructure_dlogtnorm", 0.1)

    print("done.")
    # / Read ini file

    # $ Read swap file names
    try:
        with open(swapfnfile, "r") as f:
            tokens = f.read().split()
    except OSError:
        print(f'\nCan not open file "{swapfnfile}"', end="")
        return 1
    idx = 0
    _label = tokens[idx]; idx += 1
    nfiles = int(tokens[idx]); idx += 1
    swapfnames: List[str] = []
    maxsimtime: List[float] = []
    for _q in range(nfiles):
        swapfnames.append(tokens[idx]); idx += 1
        maxsimtime.append(float(tokens[idx])); idx += 1
    print("done.")
    # / Read swap file names

    # $ Read swap files
    print("Read swap files...")
    swaprecords: List[List[SwapRecord]] = [[] for _ in range(nfiles)]
    for nf in range(nfiles):
        print(swapfnames[nf])
        try:
            inpf = open(swapfnames[nf], "r")
        except OSError:
            print(f'\nCan not open file "{swapfnames[nf]}"', end="")
            continue
        nlines = 0
        nerrors = 0
        r = SwapRecord()  # reused like the single C++ `swaprecord r;`
        with inpf:
            for line in inpf:
                nlines += 1
                if nlines == 1:
                    continue  # skip header
                line = line.rstrip("\n")
                tokens_line = line.split("\t")
                nerrors += r.set(tokens_line)
                swaprecords[nf].append(r.copy())
        print(f"Records read: {nlines - 1}, errors {nerrors}")
    print("done.")
    # / Read swap files

    # Precompute per-file position arrays for vectorised geometry.
    positions: List[np.ndarray] = []
    for nf in range(nfiles):
        if swaprecords[nf]:
            positions.append(np.array([rec.bondCenter for rec in swaprecords[nf]], dtype=float))
        else:
            positions.append(np.zeros((0, 3), dtype=float))

    # =====================================================================
    # This block combines bond snaps into subavalanches (cascades)
    # =====================================================================
    # $ Detect subavalanches
    print("Detect subavalanches...")
    swap_subavalancheID: List[List[int]] = [[] for _ in range(nfiles)]
    subavalanches: List[List[Subavalanche]] = [[] for _ in range(nfiles)]
    subavalanches_nsimultaneous: List[List[int]] = [[] for _ in range(nfiles)]

    for nf in range(nfiles):
        nsubavalanches = 0
        nswaps = len(swaprecords[nf])
        swap_subavalancheID[nf] = [-1] * nswaps

        current_subavalanches: List[List[int]] = []  # swap IDs in currently open subavalanches
        for ns in range(nswaps):
            ns_time = swaprecords[nf][ns].time
            ns_center = swaprecords[nf][ns].bondCenter

            # Calculate dt and dr for each of the current subavalanches
            nsubav = len(current_subavalanches)
            dt_last = [0.0] * nsubav
            dr_threshold_min = [-1.0] * nsubav
            dt_threshold_min = [0.0] * nsubav
            for q in range(nsubav):
                saSize = len(current_subavalanches[q])
                dt_last[q] = ns_time - swaprecords[nf][current_subavalanches[q][saSize - 1]].time
                dr_threshold_min[q] = -1.00
                for w in range(saSize):
                    other = swaprecords[nf][current_subavalanches[q][saSize - 1 - w]]
                    dt = ns_time - other.time
                    if dt > subavalanche_dt_threshold:
                        break
                    dr = float(np.linalg.norm(ns_center - other.bondCenter))
                    if dr_threshold_min[q] >= 0.00:
                        if dr_threshold_min[q] > dr:
                            dr_threshold_min[q] = dr
                            dt_threshold_min[q] = dt
                    else:
                        dr_threshold_min[q] = dr
                        dt_threshold_min[q] = dt

            # Determine the closest subavalanche
            closestSubav = -1
            dr_min = 0.00
            for q in range(nsubav):
                if dr_threshold_min[q] >= 0.00 and dr_threshold_min[q] < subavalanche_dr_threshold:
                    if closestSubav < 0:
                        dr_min = dr_threshold_min[q]
                        closestSubav = q
                    elif dr_min > dr_threshold_min[q]:
                        dr_min = dr_threshold_min[q]
                        closestSubav = q

            # Check if the new event breaks any of the subavalanche sequences
            isSubavStopped = [False] * nsubav
            for q in range(nsubav):
                if dt_last[q] > subavalanche_dt_threshold:
                    isSubavStopped[q] = True

            # Add the event to the closest subavalanche (if exists) or to a new one
            if closestSubav >= 0:
                current_subavalanches[closestSubav].append(ns)
            else:
                current_subavalanches.append([ns])
                isSubavStopped.append(False)

            # Store & forget stopped current subavalanches
            tryErasing = True
            while tryErasing:
                tryErasing = False
                for q in range(len(current_subavalanches)):
                    if isSubavStopped[q]:
                        newsa = Subavalanche()
                        newsa.subavalancheID = nsubavalanches
                        for w in range(len(current_subavalanches[q])):
                            swapID = current_subavalanches[q][w]
                            newsa.swaprecords_time_ID.append((swaprecords[nf][swapID].time, swapID))
                            swap_subavalancheID[nf][swapID] = nsubavalanches
                        subavalanches[nf].append(newsa)
                        nsubavalanches += 1
                        del current_subavalanches[q]
                        del isSubavStopped[q]
                        tryErasing = True
                        break

        # Store all lasting current subavalanches
        for q in range(len(current_subavalanches)):
            newsa = Subavalanche()
            newsa.subavalancheID = nsubavalanches
            for w in range(len(current_subavalanches[q])):
                swapID = current_subavalanches[q][w]
                newsa.swaprecords_time_ID.append((swaprecords[nf][swapID].time, swapID))
                swap_subavalancheID[nf][swapID] = nsubavalanches
            subavalanches[nf].append(newsa)
            nsubavalanches += 1

        # Calculate properties of subavalanches detected
        for nsa in range(len(subavalanches[nf])):
            sa = subavalanches[nf][nsa]
            sa.swaprecords_time_ID.sort()
            sa.nswaps = len(sa.swaprecords_time_ID)
            sa.subavTimeStart = sa.swaprecords_time_ID[0][0]
            sa.subavTimeEnd = sa.swaprecords_time_ID[sa.nswaps - 1][0]
            sa.subavDuration = sa.subavTimeEnd - sa.subavTimeStart
            nswap_P50 = sa.nswaps // 2
            if sa.nswaps % 2 == 1:
                sa.subavTimeP50 = sa.swaprecords_time_ID[nswap_P50][0]
            else:
                sa.subavTimeP50 = (sa.swaprecords_time_ID[nswap_P50][0] + sa.swaprecords_time_ID[nswap_P50 - 1][0]) * 0.50

            ids = [i for (_, i) in sa.swaprecords_time_ID]
            pts = positions[nf][ids]
            sa.subavDiameter = _max_pairwise_distance(pts) + avgBondLength
            center = pts.mean(axis=0)
            sa.subavCenterX = float(center[0])
            sa.subavCenterY = float(center[1])
            dr = center - pts
            sa.subavRg = math.sqrt(float(np.einsum("ij,ij->", dr, dr)) / float(sa.nswaps))
            if sa.nswaps <= 1:
                sa.subavI1 = sa.subavI2 = 0.00
            else:
                drx = center[0] - pts[:, 0]
                dry = center[1] - pts[:, 1]
                Ixx = float(np.sum(dry * dry))
                Iyy = float(np.sum(drx * drx))
                Ixy = float(-np.sum(drx * dry))
                d = math.sqrt((Ixx - Iyy) * (Ixx - Iyy) + 4.00 * Ixy * Ixy)
                sa.subavI1 = 0.5 * ((Ixx + Iyy) + d)
                sa.subavI2 = 0.5 * ((Ixx + Iyy) - d)

        # Mark simultaneous subavalanches
        nsub = len(subavalanches[nf])
        subavalanches_nsimultaneous[nf] = [0] * nsub
        for nsa in range(nsub):
            a = subavalanches[nf][nsa]
            cnt = 0
            for nsa2 in range(nsub):
                b = subavalanches[nf][nsa2]
                if ((a.subavTimeStart >= b.subavTimeStart and a.subavTimeStart <= b.subavTimeEnd) or
                        (a.subavTimeEnd >= b.subavTimeStart and a.subavTimeEnd <= b.subavTimeEnd) or
                        (b.subavTimeStart >= a.subavTimeStart and b.subavTimeStart <= a.subavTimeEnd) or
                        (b.subavTimeEnd >= a.subavTimeStart and b.subavTimeEnd <= a.subavTimeEnd)):
                    cnt += 1
            subavalanches_nsimultaneous[nf][nsa] = cnt

    print("done.")
    # / Detect subavalanches

    # =====================================================================
    # This block combines subavalanches (cascades) into thermal avalanches
    # =====================================================================
    # $ Detect avalanches
    print("Detect avalanches...")
    swaprecords_avalancheID: List[List[int]] = [[] for _ in range(nfiles)]
    subavalanche_avalancheID: List[List[int]] = [[] for _ in range(nfiles)]
    avalanches: List[List[Avalanche]] = [[] for _ in range(nfiles)]

    def _new_av_from_subavs(nf, subav_ids, av_id):
        """Build an avalanche from an ordered list of subavalanche IDs and tag members."""
        av = Avalanche()
        av.avalancheID = av_id
        for subavID in subav_ids:
            sa = subavalanches[nf][subavID]
            av.subavalanches_timeP50_ID.append((sa.subavTimeP50, subavID))
            av.subavalanches_timeStart_ID.append((sa.subavTimeStart, subavID))
            subavalanche_avalancheID[nf][subavID] = av_id
            for (t, sid) in sa.swaprecords_time_ID:
                av.swaprecords_time_ID.append((t, sid))
                swaprecords_avalancheID[nf][sid] = av_id
        return av

    if avalancheDetectionMechanism == 0:
        print("avalancheDetectionMechanism == 0")
        for nf in range(nfiles):
            swaprecords_avalancheID[nf] = [-1] * len(swaprecords[nf])
            subavalanche_avalancheID[nf] = [-1] * len(subavalanches[nf])

            subavSorted = sorted((subavalanches[nf][nsa].subavTimeP50, nsa)
                                 for nsa in range(len(subavalanches[nf])))

            # Detect interavalanche subavalanches
            interavalanche = [False] * len(subavSorted)
            for q in range(1, len(subavSorted)):
                dttP50 = (subavSorted[q][0] - subavSorted[q - 1][0]) / subavSorted[q][0]
                interavalanche[q] = dttP50 >= avDetection_dtt_threshold

            # Avalanche push_back cycle: accumulate subav IDs, split at interavalanche gaps
            current_subav_ids: List[int] = []
            for q in range(len(subavSorted)):
                if interavalanche[q]:
                    if len(current_subav_ids) > 0:
                        av = _new_av_from_subavs(nf, current_subav_ids, len(avalanches[nf]))
                        av.subavalanches_timeP50_ID.sort()
                        av.subavalanches_timeStart_ID.sort()
                        av.swaprecords_time_ID.sort()
                        compute_av_properties(av, swaprecords[nf], positions[nf])
                        avalanches[nf].append(av)
                    current_subav_ids = [subavSorted[q][1]]
                else:
                    current_subav_ids.append(subavSorted[q][1])
            # Store the last avalanche
            if len(current_subav_ids) > 0:
                av = _new_av_from_subavs(nf, current_subav_ids, len(avalanches[nf]))
                av.subavalanches_timeP50_ID.sort()
                av.subavalanches_timeStart_ID.sort()
                av.swaprecords_time_ID.sort()
                compute_av_properties(av, swaprecords[nf], positions[nf])
                avalanches[nf].append(av)

    elif avalancheDetectionMechanism == 1:
        print("avalancheDetectionMechanism == 1")
        for nf in range(nfiles):
            swaprecords_avalancheID[nf] = [-1] * len(swaprecords[nf])
            subavalanche_avalancheID[nf] = [-1] * len(subavalanches[nf])

            subavSorted = sorted((subavalanches[nf][nsa].subavTimeP50, nsa)
                                 for nsa in range(len(subavalanches[nf])))

            navalanches = 0
            current_avalanches: List[List[int]] = []
            for ns in range(len(subavSorted)):
                scan_P50 = subavSorted[ns][0]
                subavID = subavSorted[ns][1]
                scan = subavalanches[nf][subavID]
                nav = len(current_avalanches)
                dttP50_last = [0.0] * nav
                dr_threshold_min = [0.0] * nav
                threshold_min_defined = [False] * nav
                for q in range(nav):
                    avSize = len(current_avalanches[q])
                    lastSubavID = current_avalanches[q][avSize - 1]
                    dtP50_last = scan_P50 - subavalanches[nf][lastSubavID].subavTimeP50
                    dttP50_last[q] = dtP50_last / scan_P50
                    for w in range(avSize):
                        wID = current_avalanches[q][avSize - 1 - w]
                        other = subavalanches[nf][wID]
                        dtP50 = scan_P50 - other.subavTimeP50
                        dttP50 = dtP50 / scan_P50
                        if dttP50 > avDetection_dtt_threshold:
                            break
                        dx = scan.subavCenterX - other.subavCenterX
                        dy = scan.subavCenterY - other.subavCenterY
                        dr = math.sqrt(dx * dx + dy * dy)
                        if threshold_min_defined[q]:
                            if dr_threshold_min[q] > dr:
                                dr_threshold_min[q] = dr
                        else:
                            dr_threshold_min[q] = dr
                            threshold_min_defined[q] = True

                # Determine the closest avalanche
                closestAv = -1
                dr_min = 0.00
                for q in range(nav):
                    if threshold_min_defined[q] and dr_threshold_min[q] < avDetection_dr_threshold:
                        if closestAv < 0:
                            dr_min = dr_threshold_min[q]
                            closestAv = q
                        elif dr_min > dr_threshold_min[q]:
                            dr_min = dr_threshold_min[q]
                            closestAv = q

                # Check if the new subavalanche breaks any of the avalanche sequences
                isAvStopped = [False] * nav
                for q in range(nav):
                    if dttP50_last[q] > avDetection_dtt_threshold:
                        isAvStopped[q] = True

                if closestAv >= 0:
                    current_avalanches[closestAv].append(subavID)
                else:
                    current_avalanches.append([subavID])
                    isAvStopped.append(False)

                # Store & forget stopped current avalanches
                tryErasing = True
                while tryErasing:
                    tryErasing = False
                    for q in range(len(current_avalanches)):
                        if isAvStopped[q]:
                            av = _new_av_from_subavs(nf, current_avalanches[q], navalanches)
                            avalanches[nf].append(av)
                            navalanches += 1
                            del current_avalanches[q]
                            del isAvStopped[q]
                            tryErasing = True
                            break

            # Store all lasting current avalanches
            for q in range(len(current_avalanches)):
                av = _new_av_from_subavs(nf, current_avalanches[q], navalanches)
                avalanches[nf].append(av)
                navalanches += 1

            for av in avalanches[nf]:
                compute_av_properties(av, swaprecords[nf], positions[nf])

    elif avalancheDetectionMechanism in (2, 3, 4):
        print(f"avalancheDetectionMechanism == {avalancheDetectionMechanism}")
        for nf in range(nfiles):
            swaprecords_avalancheID[nf] = [-1] * len(swaprecords[nf])
            subavalanche_avalancheID[nf] = [-1] * len(subavalanches[nf])

            subavSorted = sorted((subavalanches[nf][nsa].subavTimeP50, nsa)
                                 for nsa in range(len(subavalanches[nf])))

            navalanches = 0
            current_avalanches_subavID: List[List[int]] = []
            current_avalanches_nswaps: List[int] = []
            for ns in range(len(subavSorted)):
                scanSubavID = subavSorted[ns][1]
                scan = subavalanches[nf][scanSubavID]
                ncav = len(current_avalanches_subavID)
                canjoin = [False] * ncav
                for na in range(ncav):
                    nsubav = len(current_avalanches_subavID[na])
                    q = nsubav - 1
                    while q >= 0 and not canjoin[na]:
                        prevSubavID = current_avalanches_subavID[na][q]
                        prev = subavalanches[nf][prevSubavID]
                        q -= 1
                        if scan.nswaps < avDetection_subavNswaps_threshold:
                            continue
                        dtP50 = scan.subavTimeP50 - prev.subavTimeP50
                        dttP50 = dtP50 / scan.subavTimeP50
                        if avalancheDetectionMechanism == 2:
                            if dttP50 > avDetection_dtt_threshold:
                                break
                        else:  # mechanisms 3 and 4
                            if dtP50 > avDetection_dt_minimum and dttP50 > avDetection_dtt_threshold:
                                break
                        dx = scan.subavCenterX - prev.subavCenterX
                        dy = scan.subavCenterY - prev.subavCenterY
                        dr = math.sqrt(dx * dx + dy * dy)
                        if avalancheDetectionMechanism == 4:
                            if dr > dtP50 * avDetection_velocity_threshold:
                                continue
                        else:  # mechanisms 2 and 3
                            if dr > avDetection_dr_threshold:
                                continue
                        canjoin[na] = True

                # Choose the biggest avalanche to fit
                biggestAvNum = -1
                biggestAvSize = -1
                for na in range(ncav):
                    if canjoin[na] and biggestAvSize < current_avalanches_nswaps[na]:
                        biggestAvSize = current_avalanches_nswaps[na]
                        biggestAvNum = na

                if biggestAvNum >= 0:
                    current_avalanches_subavID[biggestAvNum].append(scanSubavID)
                    current_avalanches_nswaps[biggestAvNum] += scan.nswaps
                else:
                    current_avalanches_subavID.append([scanSubavID])
                    current_avalanches_nswaps.append(scan.nswaps)

            # Store all lasting current avalanches
            for q in range(len(current_avalanches_subavID)):
                av = _new_av_from_subavs(nf, current_avalanches_subavID[q], navalanches)
                avalanches[nf].append(av)
                navalanches += 1

            for av in avalanches[nf]:
                compute_av_properties(av, swaprecords[nf], positions[nf])

    elif avalancheDetectionMechanism == 5:
        print("avalancheDetectionMechanism == 5")
        for nf in range(nfiles):
            swaprecords_avalancheID[nf] = [-1] * len(swaprecords[nf])
            subavalanche_avalancheID[nf] = [-1] * len(subavalanches[nf])

            subavSorted = sorted((subavalanches[nf][nsa].subavTimeP50, nsa)
                                 for nsa in range(len(subavalanches[nf])))

            navalanches = 0
            # list of [nswaps, [subavIDs...]]
            current_avalanches: List[List] = []
            for ns in range(len(subavSorted)):
                scanSubavID = subavSorted[ns][1]
                scan = subavalanches[nf][scanSubavID]
                ncav = len(current_avalanches)
                nojoin = True
                canjoin = [False] * ncav
                for na in range(ncav):
                    nsubav = len(current_avalanches[na][1])
                    q = nsubav - 1
                    while q >= 0 and not canjoin[na]:
                        prevSubavID = current_avalanches[na][1][q]
                        prev = subavalanches[nf][prevSubavID]
                        q -= 1
                        dtP50 = scan.subavTimeP50 - prev.subavTimeP50
                        dttP50 = dtP50 / scan.subavTimeP50
                        if dttP50 > avDetection_dtt_threshold:
                            break
                        dx = scan.subavCenterX - prev.subavCenterX
                        dy = scan.subavCenterY - prev.subavCenterY
                        dr = math.sqrt(dx * dx + dy * dy)
                        if dr > dtP50 * avDetection_velocity_threshold:
                            continue
                        canjoin[na] = True
                        nojoin = False

                if nojoin:
                    current_avalanches.append([scan.nswaps, [scanSubavID]])
                else:
                    for na in range(ncav):
                        if canjoin[na]:
                            current_avalanches[na][0] += scan.nswaps
                            current_avalanches[na][1].append(scanSubavID)

            # Decide what goes where, the biggest possible avalanche wins it all
            biggest_avalanches_subavID: List[List[int]] = []
            while len(current_avalanches) > 0:
                # Sort so the biggest avalanche is last (mirror C++ pair<int,vector> sort)
                current_avalanches.sort(key=lambda p: (p[0], p[1]))
                nbiggest = len(current_avalanches) - 1
                biggest = list(current_avalanches[nbiggest][1])
                biggest_avalanches_subavID.append(biggest)
                # Remove its subavalanches from all avalanches
                for subavID in biggest:
                    for entry in current_avalanches:
                        if subavID in entry[1]:
                            entry[1].remove(subavID)
                            entry[0] -= subavalanches[nf][subavID].nswaps
                if current_avalanches[nbiggest][0] == 0 and len(current_avalanches[nbiggest][1]) == 0:
                    print("+", end="")
                else:
                    print("!", end="")
                del current_avalanches[nbiggest]

            # Store biggest avalanches
            for subav_ids in biggest_avalanches_subavID:
                av = _new_av_from_subavs(nf, subav_ids, navalanches)
                avalanches[nf].append(av)
                navalanches += 1

            for av in avalanches[nf]:
                compute_av_properties(av, swaprecords[nf], positions[nf])

    print("done.")
    # / Detect avalanches

    # =====================================================================
    # Output: primary data files
    # =====================================================================

    # $ Write swap data  (<swapfile>__markup.txt)
    print("\nWriting swap data...", end="")
    for nf in range(nfiles):
        logfname = swapfnames[nf][:-4] + "__markup.txt"
        try:
            outf = open(logfname, "w")
        except OSError:
            print(f'\nCan not open file "{logfname}"', end="")
            return 1
        with outf:
            outf.write(
                "Time\tSwap #\tBond center X\tBond center Y\tBond center Z\tSurface area\t"
                "\tSubavalanche #\tSubavalanche size\tAvalanche #\tAvalanche size\n"
            )
            for ns in range(len(swaprecords[nf])):
                rec = swaprecords[nf][ns]
                subavID = swap_subavalancheID[nf][ns]
                avID = swaprecords_avalancheID[nf][ns]
                outf.write(
                    f"{fmt(rec.time, 4)}"
                    f"\t{ns + 1}"
                    f"\t{fmt(rec.bondCenter[0], 4)}"
                    f"\t{fmt(rec.bondCenter[1], 4)}"
                    f"\t{fmt(rec.bondCenter[2], 4)}"
                    f"\t{fmt(rec.surfaceArea, 4)}"
                    f"\t{subavID}"
                    f"\t{subavalanches[nf][subavID].nswaps}"
                    f"\t{avID}"
                    f"\t{avalanches[nf][avID].avNswaps}"
                    "\n"
                )
    print(" done.")
    # / Write swap data

    # $ Write subavalanche data  (<swapfile>__subavdata.txt)
    print("\nWriting subavalanche data...", end="")
    for nf in range(nfiles):
        logfname = swapfnames[nf][:-4] + "__subavdata.txt"
        try:
            outf = open(logfname, "w")
        except OSError:
            print(f'\nCan not open file "{logfname}"', end="")
            return 1
        with outf:
            outf.write(
                "Subavalanche #\tNumber of swaps\tNumber of swaps (swap direction corrected)"
                "\tSubavalanche P50 time\tSubavalanche start time\tSubavalanche end time\tSubavalanche duration"
                "\tSubavalanche mean rate\tSubavalanche diameter\tSubavalanche Rg\tSubavalanche I1\tSubavalanche I2"
                "\tSubavalanche mean density\tSubavalanche X position\tSubavalanche Y position"
                "\tSimultaneous subavalanches count\n"
            )
            for nsa in range(len(subavalanches[nf])):
                sa = subavalanches[nf][nsa]
                outf.write(
                    f"{nsa + 1}"
                    f"\t{sa.nswaps}"
                    f"\t--"
                    f"\t{fmt(sa.subavTimeP50, 4)}"
                    f"\t{fmt(sa.subavTimeStart, 4)}"
                    f"\t{fmt(sa.subavTimeEnd, 4)}"
                    f"\t{fmt(sa.subavDuration, 4)}"
                    f"\t--"
                    f"\t{fmt(sa.subavDiameter, 3)}"
                    f"\t{fmt(sa.subavRg, 3)}"
                    f"\t{fmt(sa.subavI1, 3)}"
                    f"\t{fmt(sa.subavI2, 3)}"
                    f"\t--"
                    f"\t{fmt(sa.subavCenterX, 3)}\t{fmt(sa.subavCenterY, 3)}"
                    f"\t{subavalanches_nsimultaneous[nf][nsa]}"
                    "\n"
                )
    print(" done.")
    # / Write subavalanche data

    av_header = (
        "File #\tAvalanche ID\tNumber of subavalanches\tNumber of swaps\tNumber of swaps (direction corrected)"
        "\tAvalanche start time\tAvalanche end time\tAvalanche avg time\tAvalanche P50 time"
        "\tAvalanche duration\tAvalanche Tg\tAvalanche Rg\tAvalanche Diameter\n"
    )

    def _av_row(nf, av) -> str:
        return (
            f"{nf + 1}\t{av.avalancheID + 1}"
            f"\t{av.avNsubavs}\t{av.avNswaps}\t{av.avNswapsDir}"
            f"\t{fmt(av.avTimeStart, 4)}\t{fmt(av.avTimeEnd, 4)}"
            f"\t{fmt(av.avTimeAvg, 4)}\t{fmt(av.avTimeP50, 4)}\t{fmt(av.avDuration, 4)}"
            f"\t{fmt(av.avTg, 4)}"
            f"\t{fmt(av.avRg, 3)}\t{fmt(av.avDiameter, 3)}"
            "\n"
        )

    # $ Write avalanche data  (<swapfile>__avdata.txt)
    print("\nWriting avalanche data...", end="")
    for nf in range(nfiles):
        logfname = swapfnames[nf][:-4] + "__avdata.txt"
        try:
            outf = open(logfname, "w")
        except OSError:
            print(f'\nCan not open file "{logfname}"', end="")
            return 1
        with outf:
            outf.write(av_header)
            for na in range(len(avalanches[nf])):
                outf.write(_av_row(nf, avalanches[nf][na]))
            outf.write("\n")
    print(" done.")
    # / Write avalanche data

    # $ Write all avalanche data  (<outfname>__avdata.txt)
    print("\nWriting all avalanche data...", end="")
    logf1name = outfname + "__avdata.txt"
    try:
        outf1 = open(logf1name, "w")
    except OSError:
        print(f'\nCan not open file "{logf1name}"', end="")
        return 1
    with outf1:
        outf1.write(av_header)
        for nf in range(nfiles):
            for na in range(len(avalanches[nf])):
                outf1.write(_av_row(nf, avalanches[nf][na]))
            outf1.write("\n")
    print(" done.")
    # / Write all avalanche data

    # $ Write avalanche structure data  (<outfname>_avstructure.txt)
    print("\nWriting avalanche structure data...", end="")
    logf11name = outfname + "_avstructure.txt"
    try:
        outf11 = open(logf11name, "w")
    except OSError:
        print(f'\nCan not open file "{logf11name}"', end="")
        return 1

    with outf11:
        # --- Compose linear time, absolute-t, absolute-nswaps data ---
        tmax = 0.0
        for nf in range(nfiles):
            for av in avalanches[nf]:
                if av.avNswaps >= avStructure_nswaps_threshold:
                    duration = av.avTimeEnd - av.avTimeStart
                    if tmax < duration:
                        tmax = duration
        ntsteps = int(tmax / avStructure_dt + 1)
        av_nswaps: List[List[int]] = []
        for nf in range(nfiles):
            for av in avalanches[nf]:
                if av.avNswaps >= avStructure_nswaps_threshold:
                    nswaps = [0] * ntsteps
                    nt = 0
                    t = 0.0
                    for ns in range(len(av.swaprecords_time_ID)):
                        time = av.swaprecords_time_ID[ns][0] - av.avTimeStart
                        while t <= time and nt < ntsteps:
                            nswaps[nt] = ns
                            nt += 1
                            t += avStructure_dt
                    for q in range(nt, ntsteps):
                        nswaps[q] = -1
                    av_nswaps.append(nswaps)

        outf11.write("Lin-lin absolute-absolute\n")
        outf11.write("Time")
        for q in range(len(av_nswaps)):
            outf11.write(f"\tAv {q + 1}")
        outf11.write("\n")
        for w in range(ntsteps):
            outf11.write(fmt(float(w) * avStructure_dt, 3))
            for q in range(len(av_nswaps)):
                outf11.write(f"\t{av_nswaps[q][w]}" if av_nswaps[q][w] >= 0 else "\t")
            outf11.write("\n")

        # --- Compose log time, absolute-t, absolute-nswaps data ---
        logtmax = 0.0
        for nf in range(nfiles):
            for av in avalanches[nf]:
                if av.avNswaps < avStructure_nswaps_threshold:
                    continue
                if av.avTimeEnd < analysis_time_min:
                    continue
                if analysis_time_max > 0.00 and av.avTimeStart > analysis_time_max:
                    continue
                logduration = math.log10(av.avTimeEnd) - math.log10(av.avTimeStart)
                if logtmax < logduration:
                    logtmax = logduration
        nlogtsteps = int(logtmax / avStructure_dlogt + 1)
        av_nswaps_logabs: List[List[int]] = []
        av_nswaps_logabs_max: List[int] = []
        for nf in range(nfiles):
            for av in avalanches[nf]:
                if av.avNswaps < avStructure_nswaps_threshold:
                    continue
                if av.avTimeEnd < analysis_time_min:
                    continue
                if analysis_time_max > 0.00 and av.avTimeStart > analysis_time_max:
                    continue
                nswaps = [0] * nlogtsteps
                nlogt = 0
                logt = 0.0
                for ns in range(len(av.swaprecords_time_ID)):
                    logtime = math.log10(av.swaprecords_time_ID[ns][0]) - math.log10(av.avTimeStart)
                    while logt <= logtime and nlogt < nlogtsteps:
                        nswaps[nlogt] = ns
                        nlogt += 1
                        logt += avStructure_dlogt
                for q in range(nlogt, nlogtsteps):
                    nswaps[q] = -1
                av_nswaps_logabs.append(nswaps)
                av_nswaps_logabs_max.append(av.avNswaps)

        outf11.write("\n\nLog-lin absolute-absolute\n")
        outf11.write("Log time")
        for q in range(len(av_nswaps_logabs)):
            outf11.write(f"\tAv {q + 1}")
        outf11.write("\n")
        for w in range(nlogtsteps):
            outf11.write(fmt(float(w) * avStructure_dlogt, 6))
            for q in range(len(av_nswaps_logabs)):
                outf11.write(f"\t{av_nswaps_logabs[q][w]}" if av_nswaps_logabs[q][w] >= 0 else "\t")
            outf11.write("\n")

        outf11.write("\n\nLog-lin absolute-normalized\n")
        outf11.write("Log time")
        for q in range(len(av_nswaps_logabs)):
            outf11.write(f"\tAv {q + 1}")
        outf11.write("\n")
        for w in range(nlogtsteps):
            # NB: precision stays at 6 for the normalized values here (mirrors C++).
            outf11.write(fmt(float(w) * avStructure_dlogt, 6))
            for q in range(len(av_nswaps_logabs)):
                if av_nswaps_logabs[q][w] >= 0:
                    outf11.write("\t" + fmt(float(av_nswaps_logabs[q][w]) / float(av_nswaps_logabs_max[q]), 6))
                else:
                    outf11.write("\t")
            outf11.write("\n")

        # --- Compose log time, normalized-t, absolute-nswaps data ---
        nlogtnormsteps = int(1.00 / avStructure_dlogtnorm + 1)
        av_nswaps_lognorm: List[List[int]] = []
        av_nswaps_lognorm_max: List[int] = []
        for nf in range(nfiles):
            for av in avalanches[nf]:
                if av.avNswaps < 2:
                    continue
                if av.avNswaps < avStructure_nswaps_threshold:
                    continue
                if av.avTimeEnd < analysis_time_min:
                    continue
                if analysis_time_max > 0.00 and av.avTimeStart > analysis_time_max:
                    continue
                nswaps = [0] * nlogtnormsteps
                nlogtn = 0
                logtn = 0.0
                logduration = math.log10(av.avTimeEnd) - math.log10(av.avTimeStart)
                for ns in range(len(av.swaprecords_time_ID)):
                    logtime = math.log10(av.swaprecords_time_ID[ns][0]) - math.log10(av.avTimeStart)
                    logtimen = logtime / logduration
                    while logtn <= logtimen and nlogtn < nlogtnormsteps:
                        nswaps[nlogtn] = ns
                        nlogtn += 1
                        logtn += avStructure_dlogtnorm
                for q in range(nlogtn, nlogtnormsteps):
                    nswaps[q] = -1
                av_nswaps_lognorm.append(nswaps)
                av_nswaps_lognorm_max.append(av.avNswaps)

        outf11.write("\n\nLog-lin normalized-absolute\n")
        outf11.write("Log time normalized")
        for q in range(len(av_nswaps_lognorm)):
            outf11.write(f"\tAv {q + 1}")
        outf11.write("\n")
        for w in range(nlogtnormsteps):
            outf11.write(fmt(float(w) * avStructure_dlogtnorm, 6))
            for q in range(len(av_nswaps_lognorm)):
                outf11.write(f"\t{av_nswaps_lognorm[q][w]}" if av_nswaps_lognorm[q][w] >= 0 else "\t")
            outf11.write("\n")

        outf11.write("\n\nLog-lin normalized-normalized\n")
        outf11.write("Log time normalized")
        for q in range(len(av_nswaps_lognorm)):
            outf11.write(f"\tAv {q + 1}")
        outf11.write("\n")
        for w in range(nlogtnormsteps):
            outf11.write(fmt(float(w) * avStructure_dlogtnorm, 6))
            for q in range(len(av_nswaps_lognorm)):
                if av_nswaps_lognorm[q][w] >= 0:
                    outf11.write("\t" + fmt(float(av_nswaps_lognorm[q][w]) / float(av_nswaps_lognorm_max[q]), 6))
                else:
                    outf11.write("\t")
            outf11.write("\n")
    print(" done.")
    # / Write avalanche structure data

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
