"""Recover net fluid force from native free-body acceleration samples.

This is for the prepared force-driven, unconstrained six-DOF cases only. The
initial placeholder row is excluded. Unsteady force is not called steady drag.
"""
import argparse
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from sph_segments import acceleration_rows,output_dirs

ROOT = Path(__file__).resolve().parent


def analyze(name):
    folder = ROOT / "sph" / name
    root = ET.parse(folder / "generated/Case.xml").getroot()
    body = root.find("casedef/floatings/floating")
    for tag in ("translationDOF", "rotationDOF"):
        if any(body.find(tag).get(axis) != "1" for axis in "xyz"):
            raise ValueError("Force recovery requires a fully free body.")
    if body.find("linearvel") is not None or root.find(".//forcepoints") is not None:
        raise ValueError("Imposed velocity or other external forces need separate accounting.")
    mass = float(body.find("massbody").get("value"))
    gravity = root.find("casedef/constantsdef/gravity")
    g = np.array([float(gravity.get(axis)) for axis in "xyz"])
    marker = body.get("mkbound")
    paths = [d/f"FloatingAce_mkbound_{marker}.csv" for d in output_dirs(folder)]
    rows = acceleration_rows(folder,marker)
    if not rows:
        raise ValueError("No nonzero-time native acceleration samples.")
    times = np.array([float(r["time [s]"]) for r in rows])
    acceleration = np.array([[float(r[f"acelin.{ax} [m/s^2]"]) for ax in "xyz"] for r in rows])
    if not np.isfinite(acceleration).all() or not np.all(np.diff(times) > 0):
        raise ValueError("Invalid acceleration samples or nonmonotonic timestamps.")
    schedule = body.find("linearforce")
    thrust = np.zeros_like(acceleration)
    if schedule is not None:
        st = np.array([float(q.get("time")) for q in schedule])
        for axis, name_axis in enumerate("xyz"):
            values = [float(q.get(name_axis)) for q in schedule]
            thrust[:, axis] = np.interp(times, st, values)
    fluid_force = mass * acceleration - thrust - mass * g
    result = {
        "case": name,
        "status": "native_free_body_force_recorded_not_a_resistance_validation",
        "mass_kg": mass, "time_s": times.tolist(),
        "acceleration_m_s2": acceleration.tolist(),
        "assumed_external_force_N": thrust.tolist(),
        "recovered_fluid_force_N": fluid_force.tolist(),
        "force_frame": "world/laboratory",
        "formula": "F_fluid = m * acceleration - F_external - m * gravity",
        "source": [str(path.relative_to(ROOT)) for path in paths],
        "last_external_force_N": thrust[-1].tolist(),
        "last_fluid_force_N": fluid_force[-1].tolist(),
        "limits": [
            "Net instantaneous fluid force includes startup, added mass and any wave forcing.",
            "Not a separated pressure/friction resistance measurement.",
            "Initial time-zero placeholder is excluded; no other external forces are permitted.",
            "Timestamps and acceleration are native solver samples, not second derivatives of rendered motion.",
        ],
    }
    (folder / "force_history.json").write_text(json.dumps(result, indent=2))
    with (folder / "force_history.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", "applied_x_N", "applied_y_N", "applied_z_N", "fluid_x_N", "fluid_y_N", "fluid_z_N"])
        writer.writerows(np.column_stack([times, thrust, fluid_force]).tolist())
    return {k: result[k] for k in ("case", "status", "last_external_force_N", "last_fluid_force_N")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case")
    print(json.dumps(analyze(parser.parse_args().case), indent=2))
