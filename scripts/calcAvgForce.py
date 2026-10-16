#!/usr/bin/env python3
"""Trajectory force-dipole analysis using the simulation input file.

Usage (matching the original analysis script):
    python3 trajectory_force_dipoles.py input 1 0

variable: 1=save total, 2=show total, 3=save passive, 4=show passive,
          5=save active, 6=show active, 7=save all three.
start frame: zero-based trajectory frame index.

All derivatives and dipole displacements are periodic in x and y.
WARNING: The isotropic area/repulsion stress allocation is provisional and
must be checked against the actual simulation implementation. The active term
uses zetaQ_self only; zetaQ_inter is not reconstructed without its exact law.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def read_input(path):
    defaults = dict(topology='test.top', trajectory_file='trajectory.dat',
                    lastconf_file='last_conf.dat', R=8., gamma=1.4, mu=120.,
                    **{'lambda': 2.}, kappa=1.5, zetaQ_self=0.5,
                    zetaQ_inter=0., print_conf_interval=2000)
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.replace('=', ' = ').split()
            if len(parts) < 3 or parts[1] != '=' or parts[0] not in defaults:
                continue
            key, value = parts[0], parts[2]
            defaults[key] = value if key in ('topology', 'trajectory_file', 'lastconf_file') else float(value)
    base = Path(path).resolve().parent
    for key in ('topology', 'trajectory_file', 'lastconf_file'):
        p = Path(defaults[key])
        defaults[key] = p if p.is_absolute() else base / p
    return defaults


def gradient(a):
    return ((np.roll(a, -1, axis=-1) - np.roll(a, 1, axis=-1)) / 2.,
            (np.roll(a, -1, axis=-2) - np.roll(a, 1, axis=-2)) / 2.)


def frames(path, n, start):
    """Stream frames: t = ..., b = Lx Ly, followed by N particle lines."""
    with open(path, encoding='utf-8') as fh:
        frame_no = -1
        t = None
        shape = None
        particles = []
        for line in fh:
            w = line.split()
            if not w:
                continue
            if w[0] == 't':
                if frame_no >= start and frame_no >= 0:
                    if len(particles) != n:
                        raise ValueError(f'Frame {frame_no}: expected {n} cells, got {len(particles)}')
                    yield frame_no, t, shape, particles
                frame_no += 1
                t = int(float(w[2]))
                shape = None
                particles = []
            elif w[0] == 'b':
                shape = (int(float(w[3])), int(float(w[2])))
            else:
                if frame_no >= start:
                    particles.append(w)
        if frame_no >= start and frame_no >= 0:
            if len(particles) != n:
                raise ValueError(f'Frame {frame_no}: expected {n} cells, got {len(particles)}')
            yield frame_no, t, shape, particles


def reconstruct(shape, particles):
    if shape is None:
        raise ValueError('Missing b = Lx Ly header')
    ly, lx = shape
    n = len(particles)
    phi = np.zeros((n, ly, lx), dtype=np.float64)
    com = np.zeros((n, 2))
    q = np.zeros((n, 2))
    for i, w in enumerate(particles):
        if len(w) < 11 or (len(w) - 11) % 2:
            raise ValueError(f'Malformed particle record {i}')
        com[i] = float(w[2]), float(w[3])
        q[i] = float(w[9]), float(w[10])
        for k in range(11, len(w), 2):
            site = int(float(w[k]))
            if site < 0 or site >= lx * ly:
                raise ValueError(f'Invalid site {site} for particle {i}')
            phi[i, site // lx, site % lx] = float(w[k + 1])
    return phi, com, q


def forces(phi, q, cfg):
    gamma, lam, mu, kappa = (cfg[k] for k in ('gamma', 'lambda', 'mu', 'kappa'))
    R = cfg['R']
    gx, gy = gradient(phi)
    area = np.sum(phi * phi, axis=(1, 2))
    grad2 = gx * gx + gy * gy
    # Explicit Cahn-Hilliard stress, Eq. (8) in the supplied PDF.
    p_ch = (4 * gamma / lam) * phi**2 * (1 - phi)**2 - gamma * lam * grad2
    # PROVISIONAL isotropic allocations, retained from the earlier script.
    # Replace with the simulation's actual per-cell stress construction.
    a0 = np.pi * R**2
    p_area = mu * (1 - area / a0)[:, None, None]**2 * phi**2 / np.maximum(area[:, None, None], 1e-30)
    p_rep = (kappa / lam) * phi**2 * (np.sum(phi**2, axis=0)[None, :, :] - phi**2)
    p = p_ch + p_area + p_rep
    sxx = p + 2 * gamma * lam * gx * gx
    sxy = 2 * gamma * lam * gx * gy
    syy = p + 2 * gamma * lam * gy * gy
    dsxx_dx, _ = gradient(sxx)
    dsxy_dx, dsxy_dy = gradient(sxy)
    _, dsyy_dy = gradient(syy)
    passive = np.stack((dsxx_dx + dsxy_dy, dsxy_dx + dsyy_dy), axis=-1)
    # Q=[[Q00,Q01],[Q01,-Q00]]; assumes zetaQ_self multiplies -Q grad(phi).
    zeta = cfg['zetaQ_self']
    active = -zeta * np.stack((q[:, 0, None, None] * gx + q[:, 1, None, None] * gy,
                              q[:, 1, None, None] * gx - q[:, 0, None, None] * gy), axis=-1)
    return passive, active, area


def dipoles(force, centers):
    n, ly, lx, _ = force.shape
    yy, xx = np.indices((ly, lx), dtype=float)
    tensors = np.zeros((n, 2, 2))
    net = force.sum(axis=(1, 2))
    for i in range(n):
        rx = (xx - centers[i, 0] + lx / 2) % lx - lx / 2
        ry = (yy - centers[i, 1] + ly / 2) % ly - ly / 2
        d = np.array([[np.sum(rx * force[i, :, :, 0]), np.sum(rx * force[i, :, :, 1])],
                      [np.sum(ry * force[i, :, :, 0]), np.sum(ry * force[i, :, :, 1])]])
        tensors[i] = 0.5 * (d + d.T)
    values, vectors = np.linalg.eigh(tensors)
    principal_index = np.argmax(np.abs(values), axis=1)
    principal_values = values[np.arange(n), principal_index]
    directions = np.array([vectors[i, :, principal_index[i]] for i in range(n)])
    return tensors, net, principal_values, directions


def plot(phi, com, strengths, directions, t, kind, path, show, R):
    n, ly, lx = phi.shape
    fig, ax = plt.subplots(figsize=(6, 6))
    x, y = np.arange(lx), np.arange(ly)
    for field in phi:
        if field.min() < 0.5 < field.max():
            ax.contour(x, y, field, levels=[0.5], colors=['#2e9e8b'], linewidths=0.9, alpha=0.65)
    max_abs = np.max(np.abs(strengths)) if n else 0
    for center, value, direction in zip(com, strengths, directions):
        if abs(value) < 1e-14 or max_abs == 0:
            continue
        length = 0.42 * R * np.sqrt(abs(value) / max_abs)
        color = 'tab:blue' if value > 0 else 'tab:red'
        # Plot copies across boundaries when the arrow crosses the periodic edge.
        for sx in (-lx, 0, lx):
            for sy in (-ly, 0, ly):
                c = center + [sx, sy]
                if not (-length <= c[0] <= lx + length and -length <= c[1] <= ly + length):
                    continue
                for sign in (-1, 1):
                    inner = c + sign * 0.12 * length * direction
                    outer = c + sign * length * direction
                    tail, tip = (inner, outer) if value > 0 else (outer, inner)
                    ax.annotate('', xy=tip, xytext=tail,
                                arrowprops=dict(arrowstyle='-|>', color=color, lw=1.4,
                                                mutation_scale=8))
    ax.plot([], [], color='tab:blue', label='Extensile')
    ax.plot([], [], color='tab:red', label='Contractile')
    ax.legend(loc='upper right', fontsize=8)
    ax.set(xlim=(0, lx), ylim=(0, ly), title=f'{kind.capitalize()} force dipoles, t={t}')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=250)
    if show:
        plt.show()
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='Simulation input parameter file')
    parser.add_argument('variable', type=int, choices=range(1, 8))
    parser.add_argument('start_frame', type=int, help='Zero-based starting frame')
    parser.add_argument('--output-dir', default='Video_forces')
    args = parser.parse_args()
    cfg = read_input(args.input)
    with open(cfg['topology'], encoding='utf-8') as fh:
        n = int(fh.readline().split()[0])
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    modes = {1: ('total',), 2: ('total',), 3: ('passive',), 4: ('passive',),
             5: ('active',), 6: ('active',), 7: ('passive', 'active', 'total')}[args.variable]
    show = args.variable in (2, 4, 6)
    if cfg['zetaQ_inter'] != 0:
        print('WARNING: zetaQ_inter is nonzero but its force law is unspecified; only zetaQ_self is used.')
    print('WARNING: area/repulsion isotropic stresses are provisional; check simulation source for exact forces.')
    count = 0
    for frame_no, t, shape, records in frames(cfg['trajectory_file'], n, args.start_frame):
        phi, com, q = reconstruct(shape, records)
        passive, active, area = forces(phi, q, cfg)
        for mode in modes:
            field = {'passive': passive, 'active': active, 'total': passive + active}[mode]
            tensor, net, strength, direction = dipoles(field, com)
            stem = out / f'frame_{frame_no:05d}_{mode}'
            np.savez_compressed(str(stem) + '.npz', t=t, force_density=field,
                                dipole_tensor=tensor, net_force=net, eigenvalue=strength,
                                direction=direction, centers=com, area=area, q=q)
            plot(phi, com, strength, direction, t, mode,
                 None if show else str(stem) + '.png', show, cfg['R'])
        print(f'frame={frame_no}, t={t}, particles={n}, modes={",".join(modes)}')
        count += 1
    print(f'Done: processed {count} frames. Output directory: {out}')


if __name__ == '__main__':
    main()

