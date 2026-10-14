import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import font_manager, rcParams


# Load font from file
#font_path = "/home/diogo/Fonts/Times_New_Roman_Normal.ttf"
font_path = "/home/diogo/Fonts/times.ttf"
italic_font_path = "/home/diogo/Fonts/timesi.ttf"
bold_font_path = "/home/diogo/Fonts/timesbd.ttf"

custom_font = FontProperties(fname=font_path)
legend_font = FontProperties(fname=font_path, size=12)
font = font_manager.FontProperties(fname=font_path)
fonti = font_manager.FontProperties(fname=italic_font_path)
fontbd = font_manager.FontProperties(fname=bold_font_path)

# Register font with a name
font_manager.fontManager.addfont(font_path)
font_manager.fontManager.addfont(italic_font_path)
font_manager.fontManager.addfont(bold_font_path)

# Set custom mathtext font to match your font
rcParams['mathtext.fontset'] = 'custom'

# Set roman (upright), italic, and bold versions (all Times New Roman if needed)
#rcParams['mathtext.rm'] = font.get_name()  # e.g. "Times New Roman"
#rcParams['mathtext.it'] = fonti.get_name()
#rcParams['mathtext.bf'] = fontbd.get_name()
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.bf'] = 'Times New Roman:bold'


# ======================================================================
# Usage
# ======================================================================
#
# python energy_over_time.py [input] [variable] [start line]
#
# variable == 1:
#     Write energy_vs_time.txt
#
# variable == 2:
#     Plot total free energy versus time
#
# variable == 3:
#     Plot total free energy and all three energy contributions
#
# Example:
#
#     python energy_over_time.py input 1 0
#
# ======================================================================


if len(sys.argv) != 4:
    print(sys.argv[0], "[input] [variable] [start line]")
    sys.exit(1)


input_file = sys.argv[1]
variable = int(float(sys.argv[2]))
start_configuration = int(float(sys.argv[3]))


# ======================================================================
# Default simulation parameters
# ======================================================================

seed = 4982
eq_steps = 0
steps = 1000000

print_conf_interval = 2000
dt = 1.0

R = 8.0

gamma = 1.4
mu = 120.0
llambda = 2.0
kappa = 1.5

friction = 3
omega = 0.4
friction_cell = 3

J0 = 0.005
J_Q = 0.1
zetaQ_self = 0.5
zetaQ_inter = 0
anchoring = 0
lambda_wall = 3

external_forces_file = "external.conf"
topology = "test.top"
conf_file = "start.conf"
trajectory_file = "trajectory.dat"
lastconf_file = "last_conf.dat"


# ======================================================================
# Read input file
# ======================================================================

ifile = open(input_file, "r")

for line in ifile:

    words = line.split()

    if len(words) == 0:
        continue

    if words[0] == "topology":
        topology = words[2]

    elif words[0] == "conf_file":
        conf_file = words[2]

    elif words[0] == "trajectory_file":
        trajectory_file = words[2]

    elif words[0] == "lastconf_file":
        lastconf_file = words[2]

    elif words[0] == "seed":
        seed = int(float(words[2]))

    elif words[0] == "equilibration_steps":
        eq_steps = int(float(words[2]))

    elif words[0] == "print_conf_interval":
        print_conf_interval = int(float(words[2]))

    elif words[0] == "dt":
        dt = float(words[2])

    elif words[0] == "R":
        R = float(words[2])

    elif words[0] == "J0":
        J0 = float(words[2])

    elif words[0] == "gamma":
        gamma = float(words[2])

    elif words[0] == "mu":
        mu = float(words[2])

    elif words[0] == "lambda":
        llambda = float(words[2])

    elif words[0] == "kappa":
        kappa = float(words[2])

    elif words[0] == "friction":
        friction = float(words[2])

    elif words[0] == "omega":
        omega = float(words[2])

    elif words[0] == "friction_cell":
        friction_cell = float(words[2])

    elif words[0] == "J_Q":
        J_Q = float(words[2])

    elif words[0] == "zetaQ_self":
        zetaQ_self = float(words[2])

    elif words[0] == "zetaQ_inter":
        zetaQ_inter = float(words[2])

    elif words[0] == "anchoring":
        anchoring = int(float(words[2]))

    elif words[0] == "lambda_wall":
        lambda_wall = float(words[2])

ifile.close()


# ======================================================================
# Read topology
# ======================================================================

tfile = open(topology, "r")

N = int(tfile.readline().split()[0])

species = []

line = tfile.readline()
lline = line.split()

for l in lline:
    species.append(int(l))

tfile.close()

numspecies = len(set(species))


# ======================================================================
# Read box dimensions from last configuration
# ======================================================================

lfile = open(lastconf_file, "r")

header = lfile.readline().split()
t_last = int(float(header[2]))

header = lfile.readline().split()

lx = int(float(header[2]))
ly = int(float(header[3]))

lfile.close()


# ======================================================================
# Print information
# ======================================================================

print("--------------------------------------------------")
print("Energy analysis")
print("--------------------------------------------------")
print("Number of phase fields =", N)
print("Box size               =", lx, "x", ly)
print()
print("R       =", R)
print("gamma   =", gamma)
print("mu      =", mu)
print("lambda  =", llambda)
print("kappa   =", kappa)
print("dt      =", dt)
print("--------------------------------------------------")


# ======================================================================
# Energy functions
# ======================================================================


def cahn_hilliard_energy(phi, gamma, llambda):
    """
    Cahn-Hilliard energy for one phase field.

    Eq. S1:

        F_CH,i = gamma/lambda * integral dx [
                     4 phi_i^2 (1 - phi_i)^2
                     + lambda^2 (grad phi_i)^2
                 ]

    The lattice spacing is assumed to be:

        dx = dy = 1

    Periodic boundary conditions are used.

    The gradient is evaluated with centered finite differences:

        dphi/dx = [phi(x+1)-phi(x-1)] / 2
        dphi/dy = [phi(y+1)-phi(y-1)] / 2
    """

    # ------------------------------------------------------------------
    # Double-well contribution
    # ------------------------------------------------------------------

    bulk = 4.0 * phi**2 * (1.0 - phi)**2


    # ------------------------------------------------------------------
    # Gradient contribution
    #
    # np.roll implements periodic boundary conditions.
    # ------------------------------------------------------------------

    dphi_dx = 0.5 * (
        np.roll(phi, -1, axis=1)
        - np.roll(phi, 1, axis=1)
    )

    dphi_dy = 0.5 * (
        np.roll(phi, -1, axis=0)
        - np.roll(phi, 1, axis=0)
    )

    grad_squared = (
        dphi_dx**2
        + dphi_dy**2
    )


    # ------------------------------------------------------------------
    # Spatial integral.
    #
    # Since dx = dy = 1, integral -> lattice sum.
    # ------------------------------------------------------------------

    integral = np.sum(
        bulk
        + llambda**2 * grad_squared
    )


    return (gamma / llambda) * integral



def area_energy(integral_phi2, mu, R):
    """
    Area constraint energy for one phase field.

    Eq. S2:

        F_area,i =
            mu [
                1 - 1/(pi R^2) integral phi_i^2 dx
            ]^2
    """

    A0 = pi * R**2

    return mu * (
        1.0
        - integral_phi2 / A0
    )**2


# ======================================================================
# Arrays storing energy versus time
# ======================================================================

time_conf = []

F_total_list = []
F_CH_list = []
F_area_list = []
F_rep_list = []


# ======================================================================
# Trajectory structure
# ======================================================================
#
# Your trajectory has:
#
#     t ...
#     b ...
#     phase field 0
#     phase field 1
#     ...
#     phase field N-1
#
# so there are:
#
#     N + 2
#
# lines per configuration.
#
# ======================================================================

start_line = (N + 2) * start_configuration


# ======================================================================
# Open trajectory
# ======================================================================

cfile = open(trajectory_file, "r")


# Skip configurations before requested starting point
for i in range(start_line):
    cfile.readline()


# ======================================================================
# Variables for current frame
# ======================================================================

current_time = 0

pt_num = 0

F_CH_frame = 0.0
F_area_frame = 0.0


# ======================================================================
# Arrays used for repulsive energy
#
# We use the identity:
#
#   sum_i sum_{j != i} phi_i^2 phi_j^2
#
#       =
#
#   (sum_i phi_i^2)^2 - sum_i phi_i^4
#
#
# Therefore we do not need to store every phase field simultaneously.
#
# ======================================================================

sum_phi2 = np.zeros(
    (ly, lx),
    dtype=np.float64
)

sum_phi4 = np.zeros(
    (ly, lx),
    dtype=np.float64
)


# ======================================================================
# Main trajectory loop
# ======================================================================

for line in cfile:

    words = line.split()

    if len(words) == 0:
        continue


    # ==================================================================
    # New configuration
    # ==================================================================

    if words[0] == "t":

        current_time = int(float(words[2]))

        pt_num = 0

        F_CH_frame = 0.0
        F_area_frame = 0.0

        sum_phi2.fill(0.0)
        sum_phi4.fill(0.0)


    # ==================================================================
    # Box dimensions
    # ==================================================================

    elif words[0] == "b":

        new_lx = int(float(words[2]))
        new_ly = int(float(words[3]))


        # --------------------------------------------------------------
        # Normally the box is constant.
        #
        # This also handles the case in which the box changes.
        # --------------------------------------------------------------

        if new_lx != lx or new_ly != ly:

            lx = new_lx
            ly = new_ly

            sum_phi2 = np.zeros(
                (ly, lx),
                dtype=np.float64
            )

            sum_phi4 = np.zeros(
                (ly, lx),
                dtype=np.float64
            )


    # ==================================================================
    # Phase-field line
    # ==================================================================

    else:

        # --------------------------------------------------------------
        # Construct phi_i(x,y)
        # --------------------------------------------------------------

        phi = np.zeros(
            (ly, lx),
            dtype=np.float64
        )


        # --------------------------------------------------------------
        # Same trajectory format as your original script.
        #
        # words[0:11] contain phase-field/cell information.
        #
        # Starting from words[11]:
        #
        #     site value site value site value ...
        #
        # --------------------------------------------------------------

        start_value = 11

        for i in range(start_value, len(words), 2):

            site = int(float(words[i]))
            value = float(words[i + 1])

            yy = int(site / lx)
            xx = int(site - yy * lx)

            phi[yy][xx] = value


        # ==================================================================
        # phi_i^2
        # ==================================================================

        phi2 = phi * phi


        # ==================================================================
        # Integral of phi_i^2
        #
        # dx = dy = 1
        # ==================================================================

        integral_phi2 = np.sum(phi2)


        # ==================================================================
        # Cahn-Hilliard energy
        # ==================================================================

        F_CH_i = cahn_hilliard_energy(
            phi,
            gamma,
            llambda
        )

        F_CH_frame += F_CH_i


        # ==================================================================
        # Area energy
        # ==================================================================

        F_area_i = area_energy(
            integral_phi2,
            mu,
            R
        )

        F_area_frame += F_area_i


        # ==================================================================
        # Accumulate quantities needed for repulsive energy
        # ==================================================================

        sum_phi2 += phi2

        sum_phi4 += phi2 * phi2


        # ==================================================================
        # Increment phase-field index
        # ==================================================================

        pt_num += 1


        # ==================================================================
        # If all N phase fields have been read, calculate frame energy
        # ==================================================================

        if pt_num == N:

            # --------------------------------------------------------------
            # Repulsive energy
            #
            # Eq. S3:
            #
            # F_rep =
            #
            #   sum_i sum_{j != i}
            #
            #       kappa/lambda
            #
            #       integral phi_i^2 phi_j^2 dx
            #
            #
            # Using:
            #
            # sum_i sum_{j != i} phi_i^2 phi_j^2
            #
            # =
            #
            # (sum_i phi_i^2)^2 - sum_i phi_i^4
            #
            # --------------------------------------------------------------

            repulsion_density = (
                sum_phi2**2
                - sum_phi4
            )


            F_rep_frame = (
                kappa / llambda
            ) * np.sum(repulsion_density)


            # --------------------------------------------------------------
            # Total free energy
            # --------------------------------------------------------------

            F_total_frame = (
                F_CH_frame
                + F_area_frame
                + F_rep_frame
            )


            # --------------------------------------------------------------
            # Physical simulation time
            # --------------------------------------------------------------

            physical_time = current_time * dt


            # --------------------------------------------------------------
            # Store results
            # --------------------------------------------------------------

            time_conf.append(physical_time)

            F_total_list.append(F_total_frame)
            F_CH_list.append(F_CH_frame)
            F_area_list.append(F_area_frame)
            F_rep_list.append(F_rep_frame)


            # --------------------------------------------------------------
            # Print progress
            # --------------------------------------------------------------

            print(
                physical_time,
                F_total_frame,
                F_CH_frame,
                F_area_frame,
                F_rep_frame
            )


cfile.close()


# ======================================================================
# Convert to NumPy arrays
# ======================================================================

time_conf = np.asarray(time_conf)

F_total_list = np.asarray(F_total_list)
F_CH_list = np.asarray(F_CH_list)
F_area_list = np.asarray(F_area_list)
F_rep_list = np.asarray(F_rep_list)


# ======================================================================
# variable == 1
#
# Write energies to file
# ======================================================================

if variable == 1:

    output_file = "energy_vs_time.txt"

    with open(output_file, "w") as f:

        print(
            "# time F_total F_CH F_area F_rep",
            file=f
        )

        for i in range(len(time_conf)):

            print(
                time_conf[i],
                F_total_list[i],
                F_CH_list[i],
                F_area_list[i],
                F_rep_list[i],
                file=f
            )

    print()
    print("Energy data written to:")
    print(output_file)


# ======================================================================
# variable == 2
#
# Plot TOTAL free energy versus time
# ======================================================================

elif variable == 2:

    # --------------------------------------------------------------
    # Global plot font settings
    # --------------------------------------------------------------

    plt.rcParams.update({
        "font.family": "Times New Roman",
        "font.size": 18,
        "axes.labelsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
        "mathtext.fontset": "stix",
    })


    # --------------------------------------------------------------
    # Create figure
    # --------------------------------------------------------------

    fig = plt.figure(
        figsize=(5.452423529, 4.089317647)
    )

    ax = plt.gca()


    # --------------------------------------------------------------
    # Total energy
    # --------------------------------------------------------------

    ax.plot(
        time_conf,
        F_total_list,
        "-",
        linewidth=1.5
    )


    # --------------------------------------------------------------
    # Axis labels
    # --------------------------------------------------------------

    ax.set_xlabel(
        "Time",
        fontsize=18,
        fontname="Times New Roman"
    )

    ax.set_ylabel(
        r"Total free energy",
        fontsize=18,
        fontname="Times New Roman"
    )


    # --------------------------------------------------------------
    # Tick labels
    # --------------------------------------------------------------

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=18
    )

    # Explicitly enforce Times New Roman on tick labels
    for label in ax.get_xticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(18)

    for label in ax.get_yticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(18)


    # --------------------------------------------------------------
    # Scientific notation
    # --------------------------------------------------------------

    ax.ticklabel_format(
        axis="x",
        style="sci",
        scilimits=(0, 0)
    )


    # Scientific-notation multiplier (e.g. x10^5)
    ax.xaxis.get_offset_text().set_fontname("Times New Roman")
    ax.xaxis.get_offset_text().set_fontsize(18)

    ax.yaxis.get_offset_text().set_fontname("Times New Roman")
    ax.yaxis.get_offset_text().set_fontsize(18)


    # --------------------------------------------------------------
    # Number of ticks
    # --------------------------------------------------------------

    ax.locator_params(
        axis="x",
        nbins=6
    )


    # --------------------------------------------------------------
    # Figure margins
    # --------------------------------------------------------------

    plt.subplots_adjust(
        left=0.235,
        bottom=0.235,
        right=0.95,
        top=0.95
    )


    # --------------------------------------------------------------
    # Save figure
    # --------------------------------------------------------------

    plt.savefig(
        "total_energy_vs_time.pdf",
        bbox_inches="tight"
    )

    plt.savefig(
        "total_energy_vs_time.png",
        dpi=300,
        bbox_inches="tight"
    )


    # --------------------------------------------------------------
    # Display figure
    # --------------------------------------------------------------

    plt.show()

    plt.close()


# ======================================================================
# variable == 3
#
# Plot all energy contributions
# ======================================================================

elif variable == 3:

    fig = plt.figure(
        figsize=(5.452423529, 4.089317647)
    )


    # ------------------------------------------------------------------
    # Total energy
    # ------------------------------------------------------------------

    plt.plot(
        time_conf,
        F_total_list,
        "-",
        linewidth=2.0,
        label=r"$\mathcal{F}$"
    )


    # ------------------------------------------------------------------
    # Cahn-Hilliard
    # ------------------------------------------------------------------

    plt.plot(
        time_conf,
        F_CH_list,
        "-",
        linewidth=1.5,
        label=r"$\mathcal{F}_{CH}$"
    )


    # ------------------------------------------------------------------
    # Area constraint
    # ------------------------------------------------------------------

    plt.plot(
        time_conf,
        F_area_list,
        "-",
        linewidth=1.5,
        label=r"$\mathcal{F}_{area}$"
    )


    # ------------------------------------------------------------------
    # Repulsion
    # ------------------------------------------------------------------

    plt.plot(
        time_conf,
        F_rep_list,
        "-",
        linewidth=1.5,
        label=r"$\mathcal{F}_{rep}$"
    )


    plt.xlabel(
        "Time",
        fontsize=18
    )

    plt.ylabel(
        "Free energy",
        fontsize=18
    )


    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)


    plt.ticklabel_format(
        axis="x",
        style="sci",
        scilimits=(0, 0)
    )


    plt.legend(
        frameon=False,
        fontsize=13
    )


    plt.subplots_adjust(
        left=0.235,
        bottom=0.235,
        right=0.95,
        top=0.95
    )


    plt.locator_params(
        axis="x",
        nbins=6
    )


    plt.show()

    plt.close()


# ======================================================================
# Unknown option
# ======================================================================

else:

    print()
    print("Unknown variable:", variable)
    print()
    print("Available options:")
    print("variable = 1 : write energies to energy_vs_time.txt")
    print("variable = 2 : plot total free energy vs time")
    print("variable = 3 : plot all free-energy contributions")


print("done")
