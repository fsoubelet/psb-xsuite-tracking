import logging

import xobjects as xo
import xpart as xp
import xtrack as xt

from simulation_parameters import parameters as p

# context = xo.ContextCpu()  # uncomment if no OpenMP for you
context = xo.ContextCpu(omp_num_threads="auto")  # comment out if no OpenMP for you

LOGGER = logging.getLogger(__name__)

# We load the exported thick Line and Twiss
LOGGER.info("Loading thick line from JSON")
line = xt.Line.from_json("psb_line_thick_activated_cavity.json")
twiss: xt.TwissTable = line.twiss(method="4d")

# And from this we get growth rates
LOGGER.info("Getting IBS growth rates from thick lattice")
rates = twiss.get_ibs_growth_rates(
    formalism="bjorken-mtingwa",
    total_beam_intensity=p["bunch_intensity"],
    nemitt_x=p["nemitt_x"],
    nemitt_y=p["nemitt_y"],
    sigma_delta=3.7e-3,  # madx.beam.sige / (madx.beam.beta ** 2)
    bunch_length=p["sigma_z"],
)
print(f"Growth rates (thick line) in [1/s]: {rates}")


# We load the exported thin Line and Twiss
LOGGER.info("Loading thin line from JSON")
line = xt.Line.from_json("psb_line_thin.json", _context=context)
line.build_tracker(_context=context)

# This is where you can customize your script
# Any funky business should happen here, before tracking
# - inserting kick elements
# - inserting monitors
# - etc


# Generate a Gaussian particle distribution (no mismatching, no missteering etc)
LOGGER.info("Generating Gaussian particle distribution")
particles = xp.generate_matched_gaussian_bunch(
    _context=context,
    num_particles=p["n_part"],
    total_intensity_particles=p["bunch_intensity"],
    nemitt_x=p["nemitt_x"],
    nemitt_y=p["nemitt_y"],
    sigma_z=p["sigma_z"],
    line=line,
)


# Perform tracking
LOGGER.info("Tracking")
line.track(particles, num_turns=p["num_turns"], turn_by_turn_monitor=True, with_progress=1)

# Then do your exports or post-processing here
# - exporting monitor data
# - plotting
# - etc
