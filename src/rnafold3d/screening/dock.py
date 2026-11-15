"""The search: where in the pocket does a compound sit, and how good is the best way it can sit there.

A pose is a translation, an orientation and one angle per rotatable bond. The search is the one docking
programs settled on because it works: from a random start, go downhill to the nearest local minimum, then
repeat — perturb, go downhill again, keep the result when it is better (a Metropolis step when it is not).
Many such runs from different starts, then the distinct minima are reported.

Three decisions worth stating:

* **The search runs on the grid, the answer is scored exactly.** Trilinear interpolation of a Gaussian with
  a 0.5 Å width on a 0.375 Å grid is good to about a unit in ten — fine for deciding which pose is better,
  not fine for the number that goes in a table. So every pose that is kept is scored again from the atoms
  (`score.score`), and that is the number reported.
* **The compound keeps its own shape.** Only the rotatable bonds turn; rings and bond lengths are the
  input conformer's. Ring conformers and tautomers are the preparation step's job, not the search's.
* **A run is reproducible.** Everything random comes from one seeded generator, so a given seed, compound
  and grid give exactly the same poses — the rest of the engine works the same way.

What this does not do: no receptor flexibility (the RNA is rigid, which is the standard approximation and
is also why a screen is run against two conformations rather than one), no water, no induced fit.
"""
from __future__ import annotations

import numpy as np

from . import chem, score as S

EXHAUSTIVENESS = {'quick': (4, 25), 'standard': (8, 60), 'thorough': (20, 150)}
STEP_CAP = {'quick': 60, 'standard': 150, 'thorough': 400}      # a big compound costs more, but not without end
CLUSTER_RMSD = 1.5          # Å, poses closer than this are the same answer


def random_rotation(rng):
    """A uniformly random rotation as a quaternion (Shoemake)."""
    u1, u2, u3 = rng.random(3)
    return np.array([np.sqrt(1 - u1) * np.sin(2 * np.pi * u2), np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
                     np.sqrt(u1) * np.sin(2 * np.pi * u3), np.sqrt(u1) * np.cos(2 * np.pi * u3)])


def quat_matrix(q):
    """Rotation matrix of a quaternion (normalised here, so the optimiser may leave it unnormalised)."""
    q = np.asarray(q, float)
    n = np.linalg.norm(q)
    if n < 1e-12:
        return np.eye(3)
    x, y, z, w = q / n
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


class Pose:
    """One way a compound sits in the pocket: its coordinates and what it scores."""

    def __init__(self, xyz, grid_score, exact=None, parameters=None):
        self.xyz = np.asarray(xyz, float)
        self.grid_score = float(grid_score)
        self.score = float(exact if exact is not None else grid_score)
        self.parameters = None if parameters is None else np.asarray(parameters, float)

    def __repr__(self):
        return f'<Pose score {self.score:.2f}>'


class Docker:
    """A compound, a grid and the moves between poses. Build one per compound, then `run()`.

    The expensive parts are worked out once here (the torsion tree, the atom pairs that have to be checked
    against each other inside the compound, the grid's view of the compound), because a screen docks a
    hundred thousand compounds and the per-pose work is what decides whether that finishes overnight.
    """

    def __init__(self, lig, grid, weights=None):
        self.lig = lig
        self.grid = grid
        self.weights = dict(S.WEIGHTS, **(weights or {}))
        self.prep = grid.prepare_ligand(lig, self.weights)   # the maps carry these weights already
        self.root, self.branches = chem.torsion_tree(lig)
        self.n_torsions = len(self.branches)
        self.centre = lig.xyz.mean(0)
        self.local = lig.xyz - self.centre
        self.internal = self._internal_pairs()

    # --------------------------------------------------------------------------------- the pose
    def coordinates(self, p):
        """The atoms of the pose `p` = [translation (3), quaternion (4), torsion angles]."""
        xyz = self.local
        if self.n_torsions:
            xyz = chem.apply_torsions(self.lig, p[7:], (self.root, self.branches), xyz)
            xyz = xyz - xyz.mean(0)
        return xyz @ quat_matrix(p[3:7]).T + p[:3]

    def _internal_pairs(self):
        """Atom pairs inside the compound that a torsion can bring together: different rigid pieces, and
        more than three bonds apart. Within a rigid piece nothing moves, so nothing has to be checked."""
        n = len(self.lig)
        if not self.n_torsions:
            return (np.zeros(0, int), np.zeros(0, int))
        nb = self.lig.neighbours
        far = np.ones((n, n), bool)
        for a in range(n):                                   # atoms within three bonds: bonded geometry
            front = {a}
            for _ in range(3):
                front = {k for x in front for k in nb[x]} | front
            far[a, list(front)] = False
        rot = {tuple(sorted(b)) for b in self.lig.rotatable}
        piece = np.full(n, -1)
        for a in range(n):                                   # the rigid pieces: what the torsions cut apart
            if piece[a] >= 0:
                continue
            piece[a] = a
            stack = [a]
            while stack:
                x = stack.pop()
                for k in nb[x]:
                    if piece[k] < 0 and tuple(sorted((x, k))) not in rot:
                        piece[k] = a
                        stack.append(k)
        same = piece[:, None] == piece[None, :]
        i, j = np.nonzero(np.triu(far & ~same, 1))
        return (i, j)

    def internal_energy(self, xyz):
        """Repulsion and contact inside the compound, so a torsion cannot fold it through itself."""
        i, j = self.internal
        if not len(i):
            return 0.0
        d = np.linalg.norm(xyz[i] - xyz[j], axis=1) - self.lig.radius[i] - self.lig.radius[j]
        return float(self.weights['repulsion'] * np.sum(np.where(d < 0, d * d, 0.0))
                     + self.weights['gauss1'] * np.sum(np.exp(-(d / 0.5) ** 2)))

    def energy(self, p):
        """The score of a pose, from the grid, with the compound's own strain added."""
        xyz = self.coordinates(p)
        return self.grid.energy(xyz, self.lig, self.prep) + self.internal_energy(xyz)

    # --------------------------------------------------------------------------------- the search
    def start(self, rng, tries=24):
        """A random pose inside the grid box — the best of a handful of throws.

        A pose thrown at random usually lands half inside the RNA, and a minimiser started there spends its
        steps climbing out of a clash instead of finding a site. Throwing two dozen and keeping the best is
        almost free (one batched lookup) and it is what makes the search find the deep pocket of a
        structure rather than the first dent it falls into.
        """
        lo = self.grid.lo + 2.0
        hi = np.maximum(self.grid.hi - 2.0, lo + 1e-6)
        cand = [np.concatenate([rng.uniform(lo, hi), random_rotation(rng),
                                rng.uniform(-np.pi, np.pi, self.n_torsions)]) for _ in range(tries)]
        xyz = np.array([self.coordinates(p) for p in cand])
        e = self.grid.energies(xyz, self.lig, self.prep) + self._internal_batch(xyz)
        return cand[int(np.argmin(e))]

    def perturb(self, p, rng, step=2.0):
        """A neighbouring pose: a shift, a turn, and a few torsions moved.

        One move in five is a jump rather than a step — a new orientation, and a shift of several Å — which
        is what gets a compound out of a minimum it has settled into on the wrong face of a helix.
        """
        q = np.array(p, float)
        jump = rng.random() < 0.2
        q[:3] += rng.normal(0, 4.0 if jump else step, 3)
        if jump:
            q[3:7] = random_rotation(rng)
        else:
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis) or 1.0
            ang = rng.normal(0, 0.6)
            q[3:7] = _quat_mul(np.concatenate([axis * np.sin(ang / 2), [np.cos(ang / 2)]]), q[3:7])
        if self.n_torsions:
            which = rng.random(self.n_torsions) < (0.8 if jump else 0.5)
            q[7:][which] += rng.normal(0, 2.0 if jump else 0.8, int(which.sum()))
        return q

    def _coordinate_batch(self, p, eps):
        """The pose and the poses one step along each parameter, for a forward-difference gradient.

        Built in the cheap way: a shift of the translation moves the finished atoms, a change of the
        quaternion only re-rotates them, and only a torsion costs the chain of rotations again.
        """
        local = self.local
        if self.n_torsions:
            local = chem.apply_torsions(self.lig, p[7:], (self.root, self.branches), self.local)
            local = local - local.mean(0)
        R = quat_matrix(p[3:7])
        base = local @ R.T + p[:3]
        out = [base]
        for t in range(3):
            c = base.copy()
            c[:, t] += eps
            out.append(c)
        for t in range(3, 7):
            q = np.array(p, float)
            q[t] += eps
            out.append(local @ quat_matrix(q[3:7]).T + p[:3])
        for t in range(7, len(p)):
            q = np.array(p, float)
            q[t] += eps
            l2 = chem.apply_torsions(self.lig, q[7:], (self.root, self.branches), self.local)
            out.append((l2 - l2.mean(0)) @ R.T + p[:3])
        return np.array(out)

    def energy_and_gradient(self, p, eps=0.05):
        """The score of a pose and its forward-difference gradient, in one pass over the grid.

        Every pose the gradient needs is looked up in the same call (`Grid.energies`), which is what makes
        the minimiser affordable in Python. The step is 0.05 Å rather than the usual tiny one on purpose:
        the grid is linear inside a cell and kinked at its faces, so a tiny step reports the slope of one
        cell and the minimiser zigzags, while a step of about a tenth of a cell averages over the kink. In
        a redocking test this one number was the difference between 3.9 Å and 0.8 Å. The analytic gradient — the grid's own slope, chained through
        the rotation and the torsions — is what the 100,000-compound target will need; this version is the
        one that can be checked against the energy itself, and it is fast enough for the benchmark.
        """
        X = self._coordinate_batch(np.asarray(p, float), eps)
        e = self.grid.energies(X, self.lig, self.prep) + self._internal_batch(X)
        return float(e[0]), (e[1:] - e[0]) / eps

    def _internal_batch(self, poses):
        """`internal_energy` for a stack of poses."""
        i, j = self.internal
        if not len(i):
            return np.zeros(len(poses))
        d = np.linalg.norm(poses[:, i] - poses[:, j], axis=2) - self.lig.radius[i] - self.lig.radius[j]
        return (self.weights['repulsion'] * np.where(d < 0, d * d, 0.0).sum(1)
                + self.weights['gauss1'] * np.exp(-(d / 0.5) ** 2).sum(1))

    def minimise(self, p, steps=40):
        """Downhill from a pose to the nearest minimum (L-BFGS-B on the pose parameters)."""
        from scipy.optimize import minimize
        res = minimize(self.energy_and_gradient, np.asarray(p, float), jac=True, method='L-BFGS-B',
                       options={'maxiter': steps})
        return res.x, float(res.fun)

    def run(self, exhaustiveness='standard', seed=0, n_poses=9, receptor=None, temperature=1.2):
        """Dock the compound: a list of `Pose`, best first.

        With `receptor`, every kept pose is scored again from the atoms rather than from the grid, which is
        the number to report (see the module text).
        """
        if isinstance(exhaustiveness, tuple):
            runs, steps = exhaustiveness
        else:
            runs, steps = EXHAUSTIVENESS.get(str(exhaustiveness), EXHAUSTIVENESS['standard'])
            # a bigger, floppier compound has more ways to sit: the search grows with it, as Vina's does
            steps = int(round(steps * (1.0 + 0.03 * max(0, len(self.lig) - 12) + 0.25 * self.n_torsions)))
            steps = min(steps, STEP_CAP.get(str(exhaustiveness), 150))
        rng = np.random.default_rng(seed)
        found = []
        for _r in range(runs):
            p, e = self.minimise(self.start(rng))
            best_p, best_e = p, e
            for _s in range(steps):
                q, eq = self.minimise(self.perturb(p, rng))
                if eq < e or rng.random() < np.exp(-(eq - e) / temperature):
                    p, e = q, eq
                if eq < best_e:
                    best_p, best_e = q, eq
            found.append((best_e, best_p))
        candidates = []
        for e, p in found:
            xyz = self.coordinates(p)
            exact = None if receptor is None else S.score(xyz, self.lig, receptor, self.weights) \
                + self.internal_energy(xyz)
            candidates.append(Pose(xyz, e, exact, p))
        poses = []
        for cand in sorted(candidates, key=lambda c: c.score):     # the exact score decides the order
            if any(rmsd(cand.xyz, other.xyz) < CLUSTER_RMSD for other in poses):
                continue
            poses.append(cand)
            if len(poses) >= n_poses:
                break
        return poses


def _quat_mul(a, b):
    x1, y1, z1, w1 = a
    x2, y2, z2, w2 = b
    return np.array([w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                     w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
                     w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2])


def rmsd(a, b, elements=None):
    """Distance between two poses of one compound, in place (no superposition).

    With `elements`, atoms of the same element may swap: the two oxygens of a carboxylate, the carbons of a
    phenyl ring. That is what a redocking benchmark has to measure — a pose flipped by the symmetry of the
    compound is the same pose. Matching by element is a lower bound of the proper symmetry-corrected value
    (it allows swaps the molecular graph would not), which is the side to err on only for a benchmark that
    reports it as such.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or not len(a):
        raise ValueError('two poses of one compound are needed')
    if elements is None:
        return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))
    from scipy.optimize import linear_sum_assignment
    total, n = 0.0, 0
    for e in set(elements):
        idx = [k for k, x in enumerate(elements) if x == e]
        d2 = np.sum((a[idx][:, None, :] - b[idx][None, :, :]) ** 2, axis=2)
        r, c = linear_sum_assignment(d2)
        total += float(d2[r, c].sum())
        n += len(idx)
    return float(np.sqrt(total / n))


def dock(lig, grid, receptor=None, exhaustiveness='standard', seed=0, n_poses=9, weights=None):
    """Dock one compound into one grid (the convenience wrapper around `Docker`)."""
    return Docker(lig, grid, weights).run(exhaustiveness, seed, n_poses, receptor)
