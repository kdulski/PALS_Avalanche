#!/usr/bin/env python3
"""Convert an old-style FitDetails text file into the JSON config.

Usage: fitdetails_to_json.py FitDetails_tumor [output.json]

Only the discrete-fit part is converted; the "Continous Fit" section of the
old file is ignored.
"""
import json
import re
import sys

KEYS = {
    "Type of the data": "format",
    "Name of the directory in ROOT file": "root_dir",
    "Name of the histogram in directory": "root_hist",
    "What is the [BinWidth]": "bin_width",
    "If type of data is times is it from": "device",
    "LastBin is the last bin": "min_last",
    "FirstBin is the first bin": "min_first",
    "Background estimated from the X number": "bg_bins",
    "Background estimated from the [right]": "bg_side",
    "Leave X bins": "bg_skip",
    "End of Fit range": "fit_end",
    "Start of Fit range": "fit_start",
    "First Bin Center": "range_min",
    "Last Bin Center": "range_max",
    "Do You want to fix the resolution": "fix_resolution",
    "For how many sigma Background": "bg_nsigma",
    "Number of iterations": "iterations",
    "How much partially fixed": "partial",
}


def number(s):
    v = float(s)
    return int(v) if v.is_integer() and "." not in s else v


def convert(lines):
    lines = [l.rstrip("\n") for l in lines]
    raw, res, life = {}, [], []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        key = next((k for p, k in KEYS.items() if line.startswith(p)), None)
        if key and key not in raw and i + 1 < len(lines):
            raw[key] = lines[i + 1].strip()
            i += 2
            continue
        if line.startswith("Sigma(double)"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("---"):
                tok = lines[i].split()
                if len(tok) >= 2:
                    res.append((float(tok[0]), float(tok[1])))
                i += 1
        elif line.split()[:2] == ["Lifetime", "Intensity(from"]:
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("---"):
                tok = lines[i].split()
                if len(tok) >= 3:
                    life.append((float(tok[0]), float(tok[1]), tok[2]))
                i += 1
        i += 1

    fmt = {"times": "times", "histogram": "histogram", "ROOT": "root"}[raw["format"]]
    cfg = {"input": {"format": fmt}}
    if fmt == "times":
        cfg["input"]["column"] = 1 if raw.get("device") == "oscilloscope" else 0
    if fmt == "root":
        d, h = raw.get("root_dir", "none"), raw.get("root_hist", "none")
        cfg["input"]["root_object"] = h if d == "none" else f"{d}/{h}"
    cfg["histogram"] = {"bin_width": number(raw["bin_width"]),
                        "range": [number(raw["range_min"]), number(raw["range_max"])]}
    cfg["estimation"] = {
        "min_counts_first_bin": number(raw["min_first"]),
        "min_counts_last_bin": number(raw["min_last"]),
        "background": {"side": raw["bg_side"], "bins": int(raw["bg_bins"]),
                       "skip_bins": int(raw["bg_skip"]), "n_sigma": number(raw["bg_nsigma"])},
        "fit_start_fraction": number(raw["fit_start"]),
        "fit_end": number(raw["fit_end"]),
    }
    fixed = raw.get("fix_resolution", "no") == "yes"
    cfg["resolution"] = [dict({"sigma": s, "intensity": w}, **({"fixed": True} if fixed else {}))
                         for s, w in res]
    cfg["lifetimes"] = [{"tau": t, "intensity": w, "type": ty} for t, w, ty in life]
    cfg["fit"] = {"iterations": int(raw["iterations"]),
                  "partial_variation_percent": number(raw["partial"]),
                  "lifetime_limits": [0.08, 142.0],
                  "sigma_limits": [0.01, 5.0],
                  "offset_window_sigma": 5,
                  "ops_lifetime_threshold": 0.7,
                  "first_pass_options": "RMW",
                  "options": "RM"}
    cfg["output"] = {"directory": "Results", "timestamp": True}
    return cfg


def dump(cfg):
    """JSON with one component per line (easier to edit by hand)."""
    text = json.dumps(cfg, indent=2)
    for key in ("resolution", "lifetimes"):
        items = ",\n    ".join(json.dumps(c) for c in cfg[key])
        start = text.index(f'"{key}": [')
        depth = 0
        for j in range(start + len(f'"{key}": '), len(text)):
            depth += text[j] == "["
            depth -= text[j] == "]"
            if depth == 0:
                break
        text = text[:start] + f'"{key}": [\n    {items}\n  ]' + text[j + 1:]
    # [low, high] pairs on one line
    text = re.sub(r"\[\s*([-\d.eE]+),\s*([-\d.eE]+)\s*\]", r"[\1, \2]", text)
    return text + "\n"


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        out = dump(convert(f.readlines()))
    if len(sys.argv) == 3:
        with open(sys.argv[2], "w") as f:
            f.write(out)
    else:
        sys.stdout.write(out)
