"""The camera fit (Round 7.1): framing reaches the edge, clip planes follow
the scene, and one IP3 pocket can be framed without clipping the rest."""

import numpy as np
import pytest

from ip3r.render.camera import Camera
from conftest import needs_structure


def _cam(aspect):
    cam = Camera()
    cam.aspect = aspect
    return cam


@pytest.mark.parametrize("aspect", [0.3, 1.0, 2.5])
def test_frame_reaches_the_edge_along_one_axis(aspect):
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(500, 3)) * [80.0, 30.0, 50.0]     # a flat cap
    cam = _cam(aspect).frame(pts, margin=1.06)
    ex, ey = cam.screen_extent(pts)
    assert max(ex, ey) == pytest.approx(1 / 1.06, rel=0.02)
    assert max(ex, ey) <= 1.0


def test_a_narrow_view_is_width_limited():
    """The failure Round 7.1 fixed: in a tall, narrow viewport the molecule
    fills the width and leaves the height empty. The fit is right; the
    viewport was wrong."""
    pts = np.array([[-100.0, -40, 0], [100, 40, 0]])
    cam = _cam(0.25).frame(pts)
    ex, ey = cam.screen_extent(pts)
    assert ex > 0.9 and ey < 0.2


def test_scene_sets_the_clip_planes_not_the_framed_points():
    everything = np.array([[-200.0, 0, 0], [200, 0, 0], [0, 0, -200], [0, 0, 200]])
    pocket = np.array([[190.0, -5, -5], [200, 5, 5]])
    tight = _cam(1.0).frame(pocket)
    wide = _cam(1.0).frame(pocket, scene=everything)
    assert tight.scene_radius < 20
    near, far = wide.clip_planes()
    assert far - near > 400                   # the whole tetramer stays drawn
    assert wide.distance == pytest.approx(tight.distance)
    assert max(wide.screen_extent(pocket)) == pytest.approx(1 / 1.06, rel=0.02)


@needs_structure("6DQN")
def test_neighbourhood_is_the_pocket_and_its_ligand():
    from ip3r.io.loader import load
    from ip3r.structure.ligand import ligand_sites, neighbourhood
    st = load("6DQN")
    site = ligand_sites(st)[0]
    atoms = neighbourhood(st, site, 15.0)
    assert set(site.atoms) <= set(atoms)
    d = np.linalg.norm(st.xyz[atoms][:, None] - st.xyz[site.atoms][None], axis=2).min(1)
    assert d.max() <= 15.0 + 1e-6
    assert 500 < len(atoms) < st.n_atoms // 20


def test_slab_clips_in_front_of_the_framed_points_and_follows_zoom():
    everything = np.array([[0.0, 0, -200], [0, 0, 200]])
    pocket = np.array([[-5.0, -5, -5], [5, 5, 5]])
    cam = _cam(1.0).frame(pocket, scene=everything, slab=True)
    near, _ = cam.clip_planes()
    assert near == pytest.approx(cam.distance - 5.0)       # pocket's front face
    cam.zoom(2.0)
    assert cam.clip_planes()[0] == pytest.approx(cam.distance - 5.0)
    assert _cam(1.0).frame(pocket, scene=everything).slab_front is None
