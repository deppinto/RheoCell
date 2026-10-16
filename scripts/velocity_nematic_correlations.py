#!/home/diogo/PythonEnv/bin/python

import sys
import os
import numpy as np


# ============================================================
# INPUT
# ============================================================

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} [input] [start_frame]")
    sys.exit(1)

input_file = sys.argv[1]
start_frame = int(sys.argv[2])


# ============================================================
# DEFAULT PARAMETERS
# ============================================================

dt = 1.0

topology = "test.top"
trajectory_file = "trajectory.dat"
lastconf_file = "last_conf.dat"


# ============================================================
# READ INPUT FILE
# ============================================================

with open(input_file, "r") as f:

    for line in f:

        words = line.split()

        if len(words) == 0:
            continue

        if words[0] == "topology":
            topology = words[2]

        elif words[0] == "trajectory_file":
            trajectory_file = words[2]

        elif words[0] == "lastconf_file":
            lastconf_file = words[2]

        elif words[0] == "dt":
            dt = float(words[2])


# ============================================================
# READ NUMBER OF CELLS
# ============================================================

with open(topology, "r") as f:
    N = int(f.readline().split()[0])

print("Number of cells:", N)


# ============================================================
# READ SYSTEM SIZE
# ============================================================

with open(lastconf_file, "r") as f:

    header = f.readline().split()
    header = f.readline().split()

    lx = int(float(header[2]))
    ly = int(float(header[3]))

print("System size:", lx, ly)


# ============================================================
# CORRELATION SETTINGS
# ============================================================

# Radial bin width in lattice units
dr = 1.0

# Maximum distance used for radial correlations
r_max = min(lx, ly) / 2.0

n_bins = int(r_max / dr)

r_values = (np.arange(n_bins) + 0.5) * dr


# ============================================================
# FUNCTIONS
# ============================================================

def minimum_image(dx, L):
    """
    Minimum-image displacement for periodic boundary conditions.
    """

    if dx > L / 2.0:
        dx -= L

    if dx < -L / 2.0:
        dx += L

    return dx


# ------------------------------------------------------------
# First zero crossing
# ------------------------------------------------------------

def first_zero_crossing(r, C):
    """
    Calculate the first zero crossing of C(r).

    Linear interpolation is used between the two bins
    surrounding zero.
    """

    valid = np.isfinite(C)

    rr = r[valid]
    cc = C[valid]

    if len(cc) < 2:
        return np.nan

    for i in range(1, len(cc)):

        if cc[i - 1] > 0.0 and cc[i] <= 0.0:

            r1 = rr[i - 1]
            r2 = rr[i]

            c1 = cc[i - 1]
            c2 = cc[i]

            return r1 - c1 * (r2 - r1) / (c2 - c1)

    return np.nan


# ------------------------------------------------------------
# FFT autocorrelation
# ------------------------------------------------------------

def autocorrelation_fft(field):
    """
    Periodic 2D autocorrelation of a scalar field.

    Returns

        < f(x) f(x+r) >_x

    on the full 2D periodic lattice.
    """

    F = np.fft.fft2(field)

    corr = np.fft.ifft2(
        F * np.conjugate(F)
    ).real

    corr /= field.size

    return corr


# ------------------------------------------------------------
# Radial average
# ------------------------------------------------------------

def radial_average(corr, lx, ly, dr, n_bins):
    """
    Radially average a periodic 2D correlation function.
    """

    radial_sum = np.zeros(n_bins)
    radial_count = np.zeros(n_bins)


    for iy in range(ly):

        # Periodic displacement in y
        dy = min(iy, ly - iy)

        for ix in range(lx):

            # Periodic displacement in x
            dx = min(ix, lx - ix)

            r = np.sqrt(
                dx * dx +
                dy * dy
            )

            index = int(r / dr)

            if index >= n_bins:
                continue

            radial_sum[index] += corr[iy, ix]
            radial_count[index] += 1


    radial_corr = np.full(n_bins, np.nan)

    valid = radial_count > 0

    radial_corr[valid] = (
        radial_sum[valid] /
        radial_count[valid]
    )

    return radial_corr


# ------------------------------------------------------------
# Velocity correlation
# ------------------------------------------------------------

def velocity_correlation(vx, vy):
    """
    Normalized velocity correlation

    Cv(r) =
        < delta_v(x) . delta_v(x+r) >
        --------------------------------
             < |delta_v|^2 >

    The spatially uniform velocity is removed first.
    """

    vx_fluc = vx - np.mean(vx)
    vy_fluc = vy - np.mean(vy)


    # ----------------------------------------
    # Correlation of each vector component
    # ----------------------------------------

    corr_x = autocorrelation_fft(vx_fluc)
    corr_y = autocorrelation_fft(vy_fluc)

    corr = corr_x + corr_y


    # ----------------------------------------
    # Normalize so C(0) = 1
    # ----------------------------------------

    norm = np.mean(
        vx_fluc**2 +
        vy_fluc**2
    )

    if norm > 0:
        corr /= norm
    else:
        corr[:] = np.nan


    return corr


# ------------------------------------------------------------
# Nematic correlation
# ------------------------------------------------------------

def nematic_correlation(Q00, Q01):
    """
    Normalized Q-tensor correlation

    CQ(r) =
       < Q00(x) Q00(x+r) + Q01(x) Q01(x+r) >
       ------------------------------------------------
                < Q00^2 + Q01^2 >

    We do NOT subtract the spatial mean Q.
    """

    corr_00 = autocorrelation_fft(Q00)
    corr_01 = autocorrelation_fft(Q01)

    corr = corr_00 + corr_01


    norm = np.mean(
        Q00**2 +
        Q01**2
    )

    if norm > 0:
        corr /= norm
    else:
        corr[:] = np.nan


    return corr


# ============================================================
# ARRAYS FOR CELL DATA
# ============================================================

CoMX = np.zeros(N)
CoMY = np.zeros(N)

CoMX_old = np.zeros(N)
CoMY_old = np.zeros(N)

velocity_x = np.zeros(N)
velocity_y = np.zeros(N)

Q00 = np.zeros(N)
Q01 = np.zeros(N)


# ============================================================
# LATTICE FIELDS
# ============================================================

velocity_x_grid = np.zeros((ly, lx))
velocity_y_grid = np.zeros((ly, lx))

Q00_grid = np.zeros((ly, lx))
Q01_grid = np.zeros((ly, lx))

sum_phi = np.zeros((ly, lx))


# ============================================================
# CORRELATION ACCUMULATORS
# ============================================================

Cv_sum = np.zeros(n_bins)
Cv_count = np.zeros(n_bins)

CQ_sum = np.zeros(n_bins)
CQ_count = np.zeros(n_bins)


# ============================================================
# TIME SERIES
# ============================================================

time_series = []

xi_v_time = []
xi_Q_time = []

vrms_time = []
Qrms_time = []

mean_speed_time = []


# ============================================================
# OTHER VARIABLES
# ============================================================

time_conf = []

frame_count = 0
pt_num = 0

have_velocity = False


# ============================================================
# PROCESS COMPLETED FRAME
# ============================================================

def process_frame(
    current_time,
    velocity_x_grid,
    velocity_y_grid,
    Q00_grid,
    Q01_grid,
    sum_phi
):

    # --------------------------------------------------------
    # Construct phase-field weighted lattice fields
    #
    # v(x) = sum phi_i(x) v_i / sum phi_i(x)
    #
    # Q(x) = sum phi_i(x) Q_i / sum phi_i(x)
    # --------------------------------------------------------

    vx = np.zeros_like(velocity_x_grid)
    vy = np.zeros_like(velocity_y_grid)

    q00 = np.zeros_like(Q00_grid)
    q01 = np.zeros_like(Q01_grid)


    mask = sum_phi > 1e-12


    vx[mask] = (
        velocity_x_grid[mask] /
        sum_phi[mask]
    )

    vy[mask] = (
        velocity_y_grid[mask] /
        sum_phi[mask]
    )


    q00[mask] = (
        Q00_grid[mask] /
        sum_phi[mask]
    )

    q01[mask] = (
        Q01_grid[mask] /
        sum_phi[mask]
    )


    # --------------------------------------------------------
    # Correlations
    # --------------------------------------------------------

    Cv_2d = velocity_correlation(
        vx,
        vy
    )

    CQ_2d = nematic_correlation(
        q00,
        q01
    )


    # --------------------------------------------------------
    # Radial averages
    # --------------------------------------------------------

    Cv = radial_average(
        Cv_2d,
        lx,
        ly,
        dr,
        n_bins
    )

    CQ = radial_average(
        CQ_2d,
        lx,
        ly,
        dr,
        n_bins
    )


    # --------------------------------------------------------
    # Correlation lengths
    # --------------------------------------------------------

    xi_v = first_zero_crossing(
        r_values,
        Cv
    )

    xi_Q = first_zero_crossing(
        r_values,
        CQ
    )


    # --------------------------------------------------------
    # Useful diagnostics
    # --------------------------------------------------------

    vx_fluc = vx - np.mean(vx)
    vy_fluc = vy - np.mean(vy)


    vrms = np.sqrt(
        np.mean(
            vx_fluc**2 +
            vy_fluc**2
        )
    )


    Qrms = np.sqrt(
        np.mean(
            q00**2 +
            q01**2
        )
    )


    speed = np.sqrt(
        vx**2 +
        vy**2
    )

    mean_speed = np.mean(speed)


    return (
        Cv,
        CQ,
        xi_v,
        xi_Q,
        vrms,
        Qrms,
        mean_speed
    )


# ============================================================
# OPEN TRAJECTORY
# ============================================================

start_line = (N + 2) * start_frame


with open(trajectory_file, "r") as cfile:

    # --------------------------------------------------------
    # Skip requested starting frames
    # --------------------------------------------------------

    for _ in range(start_line):
        cfile.readline()


    # --------------------------------------------------------
    # Read trajectory
    # --------------------------------------------------------

    for line in cfile:

        words = line.split()

        if len(words) == 0:
            continue


        # ====================================================
        # NEW FRAME
        # ====================================================

        if words[0] == "t":

            # ------------------------------------------------
            # Process the frame that has just finished
            # ------------------------------------------------

            if frame_count > 1 and have_velocity:

                current_time = time_conf[-1] * dt


                (
                    Cv,
                    CQ,
                    xi_v,
                    xi_Q,
                    vrms,
                    Qrms,
                    mean_speed
                ) = process_frame(
                    current_time,
                    velocity_x_grid,
                    velocity_y_grid,
                    Q00_grid,
                    Q01_grid,
                    sum_phi
                )


                # --------------------------------------------
                # Accumulate velocity correlation
                # --------------------------------------------

                valid = np.isfinite(Cv)

                Cv_sum[valid] += Cv[valid]
                Cv_count[valid] += 1


                # --------------------------------------------
                # Accumulate nematic correlation
                # --------------------------------------------

                valid = np.isfinite(CQ)

                CQ_sum[valid] += CQ[valid]
                CQ_count[valid] += 1


                # --------------------------------------------
                # Save time series
                # --------------------------------------------

                time_series.append(current_time)

                xi_v_time.append(xi_v)
                xi_Q_time.append(xi_Q)

                vrms_time.append(vrms)
                Qrms_time.append(Qrms)

                mean_speed_time.append(mean_speed)


                print(
                    "time:",
                    current_time,
                    "xi_v:",
                    xi_v,
                    "xi_Q:",
                    xi_Q
                )


            # ------------------------------------------------
            # Read new time
            # ------------------------------------------------

            t = int(float(words[2]))

            time_conf.append(t)


            # ------------------------------------------------
            # Save previous COM positions
            # ------------------------------------------------

            CoMX_old = CoMX.copy()
            CoMY_old = CoMY.copy()


            # ------------------------------------------------
            # Reset cell arrays
            # ------------------------------------------------

            CoMX = np.zeros(N)
            CoMY = np.zeros(N)

            velocity_x = np.zeros(N)
            velocity_y = np.zeros(N)

            Q00 = np.zeros(N)
            Q01 = np.zeros(N)


            # ------------------------------------------------
            # Reset lattice fields
            # ------------------------------------------------

            velocity_x_grid = np.zeros((ly, lx))
            velocity_y_grid = np.zeros((ly, lx))

            Q00_grid = np.zeros((ly, lx))
            Q01_grid = np.zeros((ly, lx))

            sum_phi = np.zeros((ly, lx))


            pt_num = 0
            frame_count += 1


        # ====================================================
        # BOX INFORMATION
        # ====================================================

        elif words[0] == "b":

            lx = int(float(words[2]))
            ly = int(float(words[3]))


        # ====================================================
        # CELL DATA
        # ====================================================

        else:

            # ------------------------------------------------
            # Cell COM
            # ------------------------------------------------

            CoMX[pt_num] = float(words[2])
            CoMY[pt_num] = float(words[3])


            # ------------------------------------------------
            # Cell velocity
            # ------------------------------------------------

            if frame_count > 1:

                dx = (
                    CoMX[pt_num] -
                    CoMX_old[pt_num]
                )

                dy = (
                    CoMY[pt_num] -
                    CoMY_old[pt_num]
                )


                dx = minimum_image(
                    dx,
                    lx
                )

                dy = minimum_image(
                    dy,
                    ly
                )


                delta_t = (
                    time_conf[-1] -
                    time_conf[-2]
                ) * dt


                if delta_t > 0:

                    velocity_x[pt_num] = (
                        dx / delta_t
                    )

                    velocity_y[pt_num] = (
                        dy / delta_t
                    )

                    have_velocity = True


            # ------------------------------------------------
            # Correct nematic tensor extraction
            #
            # words[9]  = Q00
            # words[10] = Q01
            # ------------------------------------------------

            Q00[pt_num] = float(words[9])
            Q01[pt_num] = float(words[10])


            # ------------------------------------------------
            # Read phase field of this cell
            #
            # Entries are:
            #
            # site_1 phi_1 site_2 phi_2 ...
            #
            # starting at words[11]
            # ------------------------------------------------

            for i in range(11, len(words), 2):

                site = int(float(words[i]))
                value = float(words[i + 1])


                yy = int(site / lx)
                xx = site - yy * lx


                # --------------------------------------------
                # Phase-field weighted velocity
                # --------------------------------------------

                velocity_x_grid[yy, xx] += (
                    value *
                    velocity_x[pt_num]
                )

                velocity_y_grid[yy, xx] += (
                    value *
                    velocity_y[pt_num]
                )


                # --------------------------------------------
                # Phase-field weighted Q tensor
                # --------------------------------------------

                Q00_grid[yy, xx] += (
                    value *
                    Q00[pt_num]
                )

                Q01_grid[yy, xx] += (
                    value *
                    Q01[pt_num]
                )


                # --------------------------------------------
                # Total phase field
                # --------------------------------------------

                sum_phi[yy, xx] += value


            pt_num += 1


# ============================================================
# PROCESS FINAL FRAME
# ============================================================

if frame_count > 1 and have_velocity:

    current_time = time_conf[-1] * dt


    (
        Cv,
        CQ,
        xi_v,
        xi_Q,
        vrms,
        Qrms,
        mean_speed
    ) = process_frame(
        current_time,
        velocity_x_grid,
        velocity_y_grid,
        Q00_grid,
        Q01_grid,
        sum_phi
    )


    valid = np.isfinite(Cv)

    Cv_sum[valid] += Cv[valid]
    Cv_count[valid] += 1


    valid = np.isfinite(CQ)

    CQ_sum[valid] += CQ[valid]
    CQ_count[valid] += 1


    time_series.append(current_time)

    xi_v_time.append(xi_v)
    xi_Q_time.append(xi_Q)

    vrms_time.append(vrms)
    Qrms_time.append(Qrms)

    mean_speed_time.append(mean_speed)


# ============================================================
# TIME-AVERAGED CORRELATIONS
# ============================================================

Cv_avg = np.full(n_bins, np.nan)
CQ_avg = np.full(n_bins, np.nan)


valid = Cv_count > 0

Cv_avg[valid] = (
    Cv_sum[valid] /
    Cv_count[valid]
)


valid = CQ_count > 0

CQ_avg[valid] = (
    CQ_sum[valid] /
    CQ_count[valid]
)


# ============================================================
# CORRELATION LENGTH OF TIME-AVERAGED CURVES
# ============================================================

xi_v_avg = first_zero_crossing(
    r_values,
    Cv_avg
)

xi_Q_avg = first_zero_crossing(
    r_values,
    CQ_avg
)


print()
print("========================================")
print("TIME-AVERAGED RESULTS")
print("========================================")
print("Velocity correlation length:", xi_v_avg)
print("Nematic correlation length :", xi_Q_avg)
print("========================================")


# ============================================================
# SAVE SPATIAL CORRELATIONS
# ============================================================

with open("spatial_correlations.txt", "w") as f:

    print(
        "# r Cv CQ",
        file=f
    )

    for i in range(n_bins):

        print(
            r_values[i],
            Cv_avg[i],
            CQ_avg[i],
            file=f
        )


# ============================================================
# SAVE CORRELATION LENGTHS
# ============================================================

with open("correlation_lengths.txt", "w") as f:

    print(
        "# xi_velocity xi_Q",
        file=f
    )

    print(
        xi_v_avg,
        xi_Q_avg,
        file=f
    )


# ============================================================
# SAVE TIME SERIES
# ============================================================

with open("correlation_time_series.txt", "w") as f:

    print(
        "# time xi_velocity xi_Q vrms Qrms mean_speed",
        file=f
    )

    for i in range(len(time_series)):

        print(
            time_series[i],
            xi_v_time[i],
            xi_Q_time[i],
            vrms_time[i],
            Qrms_time[i],
            mean_speed_time[i],
            file=f
        )


print("done")
