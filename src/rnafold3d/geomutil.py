"""Small, well-tested geometry primitives shared by the v2 engine."""
from __future__ import annotations

import math

import numpy as np


def angle(a, b, c):
    """Bond angle a-b-c in degrees."""
    u = np.asarray(a, float) - b
    v = np.asarray(c, float) - b
    cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def dihedral(a, b, c, d):
    """IUPAC dihedral a-b-c-d in degrees, range (-180, 180]."""
    b0 = np.asarray(a, float) - b
    b1 = np.asarray(c, float) - b
    b2 = np.asarray(d, float) - c
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def kabsch(src, dst):
    """Return (R, t) minimizing |src @ R.T + t - dst| (proper rotation)."""
    src = np.asarray(src, float)
    dst = np.asarray(dst, float)
    cs = src.mean(axis=0)
    cd = dst.mean(axis=0)
    h = (src - cs).T @ (dst - cd)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, cd - cs @ r.T


def rmsd_after_fit(a, b):
    r, t = kabsch(a, b)
    return float(np.sqrt(np.mean(np.sum((np.asarray(a) @ r.T + t - b) ** 2, axis=1))))


def place(a, b, c, bond, ang_deg, tors_deg):
    """NeRF: position of atom d given a, b, c, |cd|, angle(b,c,d) and dihedral(a,b,c,d).

    Scalar implementation (about 30x faster than vectorised numpy for single atoms)."""
    ax, ay, az = a.tolist() if hasattr(a, 'tolist') else a
    bx, by, bz = b.tolist() if hasattr(b, 'tolist') else b
    cx, cy, cz = c.tolist() if hasattr(c, 'tolist') else c
    ux, uy, uz = cx - bx, cy - by, cz - bz
    lu = math.sqrt(ux * ux + uy * uy + uz * uz)
    ux, uy, uz = ux / lu, uy / lu, uz / lu
    vx, vy, vz = bx - ax, by - ay, bz - az
    nx, ny, nz = vy * uz - vz * uy, vz * ux - vx * uz, vx * uy - vy * ux
    ln = math.sqrt(nx * nx + ny * ny + nz * nz)
    nx, ny, nz = nx / ln, ny / ln, nz / ln
    mx, my, mz = ny * uz - nz * uy, nz * ux - nx * uz, nx * uy - ny * ux
    th = math.radians(ang_deg)
    ph = math.radians(tors_deg)
    d0 = -bond * math.cos(th)
    st = bond * math.sin(th)
    d1 = st * math.cos(ph)
    d2 = st * math.sin(ph)
    return np.array([cx + d0 * ux + d1 * mx + d2 * nx,
                     cy + d0 * uy + d1 * my + d2 * ny,
                     cz + d0 * uz + d1 * mz + d2 * nz])


def frame(p0, p1, p2):
    """Orthonormal frame (rows e1, e2, e3) with origin p1: e1 along p0-p1, e2 in the plane towards p2."""
    x0, y0, z0 = p0.tolist() if hasattr(p0, 'tolist') else p0
    x1, y1, z1 = p1.tolist() if hasattr(p1, 'tolist') else p1
    x2, y2, z2 = p2.tolist() if hasattr(p2, 'tolist') else p2
    ex, ey, ez = x0 - x1, y0 - y1, z0 - z1
    le = math.sqrt(ex * ex + ey * ey + ez * ez)
    ex, ey, ez = ex / le, ey / le, ez / le
    wx, wy, wz = x2 - x1, y2 - y1, z2 - z1
    d = wx * ex + wy * ey + wz * ez
    wx, wy, wz = wx - d * ex, wy - d * ey, wz - d * ez
    lw = math.sqrt(wx * wx + wy * wy + wz * wz)
    wx, wy, wz = wx / lw, wy / lw, wz / lw
    return np.array([[ex, ey, ez], [wx, wy, wz],
                     [ey * wz - ez * wy, ez * wx - ex * wz, ex * wy - ey * wx]])


def to_frame(points, p0, p1, p2):
    """Express points in the frame(p0, p1, p2) (origin p1)."""
    F = frame(p0, p1, p2)
    return (np.asarray(points, float) - p1) @ F.T


def from_frame(local, p0, p1, p2):
    F = frame(p0, p1, p2)
    return local @ F + np.asarray(p1, float)


def rotation_about_axis(axis, angle_rad):
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    q = 1 - c
    return np.array([[c + x * x * q, x * y * q - z * s, x * z * q + y * s],
                     [y * x * q + z * s, c + y * y * q, y * z * q - x * s],
                     [z * x * q - y * s, z * y * q + x * s, c + z * z * q]])


def pseudorotation(ring):
    """Altona-Sundaralingam phase and amplitude (deg) from (C4', O4', C1', C2', C3') coordinates."""
    c4, o4, c1, c2, c3 = ring
    nu0 = dihedral(c4, o4, c1, c2)
    nu1 = dihedral(o4, c1, c2, c3)
    nu2 = dihedral(c1, c2, c3, c4)
    nu3 = dihedral(c2, c3, c4, o4)
    nu4 = dihedral(c3, c4, o4, c1)
    num = (nu4 + nu1) - (nu3 + nu0)
    den = 2.0 * nu2 * (math.sin(math.radians(36)) + math.sin(math.radians(72)))
    phase = math.degrees(math.atan2(num, den)) % 360
    amp = nu2 / math.cos(math.radians(phase)) if abs(math.cos(math.radians(phase))) > 1e-6 else float('nan')
    return phase, amp


def circular_mean(deg):
    rad = np.radians(np.asarray(deg, float))
    return float(np.degrees(np.arctan2(np.sin(rad).mean(), np.cos(rad).mean())))


def wrap(deg):
    return (np.asarray(deg) + 180.0) % 360.0 - 180.0
