import numpy as np

from structcooker.design.ss import assign8, max_run, to_three


def build(phi, psi, n):
    """Backbone N CA C O of an n-residue chain with constant phi / psi (NeRF)."""
    ang = {"N-CA": 1.458, "CA-C": 1.525, "C-N": 1.329, "CA-C-N": 116.2, "C-N-CA": 121.7, "N-CA-C": 111.0}
    pts = [np.array([0.0, 0.0, 0.0]), np.array([1.458, 0.0, 0.0])]
    a = np.radians(111.0)
    pts.append(pts[1] + 1.525 * np.array([-np.cos(a), np.sin(a), 0.0]))
    seq = [("CA-C-N", "C-N", psi), ("C-N-CA", "N-CA", 180.0), ("N-CA-C", "CA-C", phi)]
    k = 0
    while len(pts) < 3 * n:
        ba, bl, t = seq[k % 3]
        a_, b_, c_ = pts[-3], pts[-2], pts[-1]
        bc = (c_ - b_) / np.linalg.norm(c_ - b_)
        nn = np.cross(b_ - a_, bc)
        nn /= np.linalg.norm(nn)
        m = np.stack([bc, np.cross(nn, bc), nn], axis=1)
        th, ph = np.radians(180 - ang[ba]), np.radians(t)
        d = np.array([np.cos(th), np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)])
        pts.append(c_ + ang[bl] * (m @ d))
        k += 1
    p = np.array(pts[: 3 * n]).reshape(n, 3, 3)
    o = np.zeros((n, 3))
    for i in range(n):
        v = p[i, 2] - p[i, 1]
        u = p[i, 1] - p[i, 0]
        w = v / np.linalg.norm(v) - u / np.linalg.norm(u)
        o[i] = p[i, 2] + 1.23 * w / np.linalg.norm(w) if i == n - 1 else p[i, 2] + 1.23 * _perp(p, i)
    return np.concatenate([p, o[:, None]], axis=1)


def _perp(p, i):
    c_o = (p[i, 2] - p[i, 1]) / np.linalg.norm(p[i, 2] - p[i, 1])
    nxt = p[i + 1, 0] - p[i, 2]
    nxt /= np.linalg.norm(nxt)
    v = c_o - nxt
    return v / np.linalg.norm(v)


def run(x):
    n = len(x)
    return assign8(x, np.zeros(n, dtype=int), np.arange(n), np.zeros(n, dtype=bool), np.ones(n, dtype=bool))


def test_alpha_helix():
    s = run(build(-57.0, -47.0, 30))
    assert (s[3:-3] == "H").all()


def test_extended_strand_is_not_helix():
    s = to_three(run(build(-120.0, 130.0, 20)))
    assert "H" not in s


def test_max_run():
    assert max_run(np.array([1, 1, 0, 1, 1, 1, 0], dtype=bool)) == 3
