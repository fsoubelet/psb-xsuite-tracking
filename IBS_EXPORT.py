import logging

import tfs
import xpart as xp
import xtrack as xt
from cpymad.madx import Madx

from simulation_parameters import parameters as p

LOGGER = logging.getLogger(__name__)


def get_sire_twiss(line: xt.Line) -> tfs.TfsDataFrame:
    """
    Get a SIRE-compatible `TfsDataFrame` to write to disk from the line or
    a twiss (4d is fine) of the line.

    Parameters
    ----------
    line : xt.Line, optional
        The line itself, it is expected a reference particle has been set.

    Returns
    -------
    df : tfs.TfsDataFrame
        A TfsDataframe ready to write to disk, that can be read by SIRE.
        Write it with: tfs.write("path/to/location.tfs", df)
    """
    # Get a Twiss and convert to pandas.DataFrame
    twiss = line.twiss(method="4d")
    df = twiss.to_pandas()

    # Start creating our own
    result = tfs.TfsDataFrame()
    result["NAME"] = df.name.str.upper()
    result["S"] = df.s
    # Now we need to get elements lengths which is not always straightforward in Xsuite
    # Best is to get S coordinate diff, shifted by -1 (so first element has a length -> last one has NaN).
    # This is ok because Xsuite inserts an "_end_point" marker at end of line, which has 0 length
    result["L"] = result.S.diff().shift(-1).fillna(0)
    result["BETX"] = df.betx
    result["ALFX"] = df.alfx
    result["MUX"] = df.mux
    result["BETY"] = df.bety
    result["ALFY"] = df.alfy
    result["MUY"] = df.muy
    result["DX"] = df.dx
    result["DPX"] = df.dpx
    result["DY"] = df.dy
    result["DPY"] = df.dpy

    # Now it needs headers, not sure if SIRE actually cares
    result.headers = {
        "TITLE": "TWISS FOR SIRE",
        "MASS [GeV]": line.particle_ref.mass0 / 1e9,
        "CHARGE": line.particle_ref.charge[0],
        "ENERGY": line.particle_ref.energy0[0] / 1e9,
        "PC": line.particle_ref.p0c[0] / 1e9,
        "GAMMA": line.particle_ref.gamma0[0],
        "KBUNCH": 1.0,
        "LENGTH": twiss.circumference,
        "Q1": twiss.qx,
        "Q2": twiss.qy,
        "DQ1": twiss.dqx,
        "DQ2": twiss.dqy,
    }

    return result


# This first parts gets the lattice directly through MAD-X from the
# scripts in the psb folder (can be updated)
LOGGER.info("Getting the lattice from MAD-X")
madx = Madx(command_log="commands.madx", stdout=False)
madx.globals.QH = p["qx_ini"]
madx.globals.QV = p["qy_ini"]
madx.chdir("psb")  # scripts in here call direct file names so we need to be there
madx.call("psb_flat_bottom.madx")  # a single cavity is installed close to the main one (br.c02)

# Now we load this as a Line in xtrack
LOGGER.info("Loading the lattice in xtrack")
line = xt.Line.from_madx_sequence(
    madx.sequence.psb1,
    deferred_expressions=True,
    install_apertures=True,
    enable_align_errors=True,
    allow_thick=True,
)
line.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, gamma0=madx.sequence.psb1.beam.gamma)
madx.exit()

# And we configure the bend model to be exact (recommended for short bends, appropriate for small rings)
LOGGER.debug("Configuring the bend model to full")
line.configure_bend_model(core="full", edge="full")

# Get a twiss, export to json
LOGGER.info("Getting twiss and exporting to json")
twiss = line.twiss(method="4d")
line.to_json("psb_line_thick.json")

# Now to configure for longitudinal motion (single RF, no acceleration)
LOGGER.info("Configuring Single RF for longitudinal motion")
line.element_refs["br.c02"].voltage = 0.008 * 1e6  # voltage 8kV

# Get a twiss, export to json
LOGGER.info("Getting twiss and exporting to json")
twiss2 = line.twiss(method="6d")
line.to_json("psb_line_thick_activated_cavity.json")


# Getting SIRE-compatible twiss
LOGGER.info("Getting SIRE-compatible twiss and exporting")
siretwiss = get_sire_twiss(line)
location = "for_sire.tfs"
tfs.write(location, siretwiss)

# Potentially we can slice the line (not strictly needed, but makes tracking faster)
# Deactivate painting bump, chicane and correction
LOGGER.info("Slicing the line")
line.vars["on_painting_bump"] = 0
line.vars["on_chicane_k0"] = 0
line.vars["on_chicane_k2"] = 0
line.vars["on_chicane_tune_corr"] = 0
line.vars["on_chicane_beta_corr"] = 0

slices = p["slices"]
line.slice_thick_elements(
    slicing_strategies=[
        xt.Strategy(slicing=xt.Teapot(slices)),
        xt.Strategy(slicing=xt.Teapot(slices), element_type=xt.Bend),
        xt.Strategy(slicing=xt.Teapot(slices), element_type=xt.Quadrupole),
        xt.Strategy(slicing=xt.Teapot(slices), name=r"bi.*bsw.*"),
    ]
)

# Let's rematch the tunes as they can drift a bit with slicing
LOGGER.info("Rematching the tunes")
qx_target = p["qx_ini"]
qy_target = p["qy_ini"]
line.match(
    vary=[xt.Vary("kbrqf", step=1e-8), xt.Vary("kbrqd", step=1e-8)],
    targets=[xt.Target("qx", qx_target, tol=1e-5), xt.Target("qy", qy_target, tol=1e-5)],
)

# Reactivate painting bump, chicane and correction
line.vars["on_chicane_k0"] = p["on_chicane_k0"]
line.vars["on_chicane_k2"] = p["on_chicane_k2"]
line.vars["on_chicane_tune_corr"] = p["on_chicane_tune_corr"]
line.vars["on_chicane_beta_corr"] = p["on_chicane_beta_corr"]
line.vars["on_painting_bump"] = p["on_painting_bump"]

# Get a twiss, export to json
LOGGER.info("Getting twiss and exporting to json")
twiss2 = line.twiss()
line.to_json("psb_line_thin.json")

# We are done
print("Ran MAD-X scripts at 'commands.madx'")
print(f"SIRE-compatible TWISS at: '{location}'")
print("Line versions saved to JSON:")
print(" - Initial (thick) line from MAD-X: 'psb_line_thick.json'")
print(" - Active RF (thick) line: 'psb_line_thick_activated_cavity.json'")
print(" - Thin line: 'psb_line_thin.json')")
