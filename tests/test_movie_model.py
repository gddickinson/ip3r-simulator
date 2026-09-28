"""Movies of the 3-D view: the frame plan loops without a repeated frame,
the GIF keeps every frame and its colours, the MP4 plays at the size the
codec allows, and a QImage's row padding never leaks into a frame."""

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.ui.movie_model import (MovieSpec, delay_ms, fit_size, formats_available,
                                 gif_delay_ms, pingpong, plan, write_gif, write_mp4,
                                 write_movie, write_webp)


def _spec(source, **kw):
    return MovieSpec(source=source, frames=kw.pop("frames", 12), fps=10.0, **kw)


def test_pingpong_by_hand():
    assert pingpong(5) == [0, 1, 2, 3, 4, 3, 2, 1]          # back into 0, no repeat
    assert pingpong(3, hold=2) == [0, 0, 0, 1, 2, 2, 2, 1]
    with pytest.raises(ValueError):
        pingpong(1)


def test_turntable_closes_one_step_short():
    frames = plan(_spec("turntable", frames=8))
    angles = [f.angle for f in frames]
    assert angles == [45.0 * k for k in range(8)]           # 360 would repeat 0
    assert all(f.state is None for f in frames)


def test_mode_cycle_and_optional_turn():
    still = plan(_spec("mode", frames=4))
    assert np.allclose([f.state for f in still], [0, np.pi / 2, np.pi, 3 * np.pi / 2])
    assert all(f.angle == 0.0 for f in still)
    turning = plan(_spec("mode", frames=4, turn=True))
    assert [f.angle for f in turning] == [0.0, 90.0, 180.0, 270.0]


def test_transition_uses_the_solved_frames():
    frames = plan(_spec("transition", hold=1), n_solved=4)
    assert [f.state for f in frames] == [0, 0, 1, 2, 3, 3, 2, 1]
    with pytest.raises(ValueError):
        plan(_spec("transition"))                           # frame count unknown


def test_spec_refuses_nonsense():
    for bad in (dict(source="spin"), dict(source="mode", fmt="avi"),
                dict(source="mode", frames=1), dict(source="mode", hold=-1)):
        with pytest.raises(ValueError):
            MovieSpec(**{"frames": 12, "fps": 10.0, **bad})


def test_fit_size_never_enlarges_and_evens():
    assert fit_size(1000, 500, 720) == (720, 360)
    assert fit_size(300, 201, 720) == (300, 201)            # smaller: untouched
    assert fit_size(301, 201, 720, even=True) == (300, 200)


def test_defaults_make_a_valid_movie():
    spec = MovieSpec("turntable", int(_P.value("movie.turntable_frames")),
                     float(_P.value("movie.fps")), hold=int(_P.value("movie.hold_frames")),
                     max_width=int(_P.value("movie.max_width")))
    assert len(plan(spec)) == _P.value("movie.turntable_frames")
    assert 1.0 <= 360.0 / len(plan(spec)) <= 10.0           # a smooth turn


def _frames(n=6, h=20, w=30):
    """Frames of distinct flat colours, so each can be identified after encoding."""
    rng = np.random.default_rng(0)
    colours = rng.integers(0, 256, size=(n, 3), dtype=np.uint8)
    return [np.broadcast_to(c, (h, w, 3)).copy() for c in colours], colours


def test_gif_keeps_every_frame_and_colour(tmp_path):
    from PIL import Image
    frames, colours = _frames()
    path = write_gif(tmp_path / "m.gif", frames, fps=10.0)
    with Image.open(path) as im:
        assert im.n_frames == len(frames)
        assert im.info["loop"] == 0 and im.info["duration"] == 100
        for k, c in enumerate(colours):
            im.seek(k)
            got = np.asarray(im.convert("RGB"))[10, 15]
            assert np.abs(got.astype(int) - c).max() <= 8, (k, got, c)


def test_gif_timing_is_what_it_can_hold():
    assert gif_delay_ms(20) == 50 and gif_delay_ms(10) == 100
    assert gif_delay_ms(15) == 70          # 66.7 ms is not a whole centisecond
    assert gif_delay_ms(1000) == 10        # the shortest delay a GIF holds


def test_held_frames_keep_their_time(tmp_path):
    """A transition's held ends are identical frames: stored once, their
    delays added, so the movie plays for exactly its plan."""
    from PIL import Image
    frames, _ = _frames(n=3)
    held = [frames[0]] * 4 + frames[1:]
    path = write_gif(tmp_path / "h.gif", held, fps=20.0)
    with Image.open(path) as im:
        durations = []
        for k in range(im.n_frames):
            im.seek(k)
            durations.append(im.info["duration"])
    assert sum(durations) == len(held) * 50
    assert durations[0] == 4 * 50


def test_webp_keeps_every_frame_colour_and_time(tmp_path):
    from PIL import Image
    if formats_available()["webp"]:
        pytest.skip("Pillow without WebP")
    frames, colours = _frames(n=4)
    path = write_webp(tmp_path / "m.webp", frames, fps=15.0, quality=90)
    assert delay_ms(15.0, "webp") == 67                 # WebP holds milliseconds
    with Image.open(path) as im:
        assert im.n_frames == 4
        total = 0
        for k, c in enumerate(colours):
            im.seek(k)
            got = np.asarray(im.convert("RGB"))[10, 15]      # loads: then the delay
            total += im.info["duration"]
            assert np.abs(got.astype(int) - c).max() <= 12, (k, got, c)
    assert total == 4 * 67


def test_write_movie_webp_uses_the_registered_quality(tmp_path):
    if formats_available()["webp"]:
        pytest.skip("Pillow without WebP")
    frames, _ = _frames(n=2)
    out = write_movie(tmp_path / "clip", frames, _spec("mode", fmt="webp"))
    assert out.suffix == ".webp" and out.stat().st_size > 0


def test_mp4_frames_and_even_size(tmp_path):
    if formats_available()["mp4"]:
        pytest.skip("no imageio-ffmpeg")
    import imageio_ffmpeg
    frames, _ = _frames(n=5, h=21, w=31)                    # odd: cropped to 30 x 20
    path = write_mp4(tmp_path / "m.mp4", frames, fps=10.0)
    reader = imageio_ffmpeg.read_frames(str(path))
    meta = next(reader)
    assert tuple(meta["size"]) == (30, 20)
    assert sum(1 for _ in reader) == 5


def test_write_movie_adds_the_suffix(tmp_path):
    frames, _ = _frames(n=2)
    out = write_movie(tmp_path / "clip", frames, _spec("mode", fmt="gif"))
    assert out.suffix == ".gif" and out.exists()


def test_image_array_drops_row_padding():
    QtGui = pytest.importorskip("PyQt6.QtGui")
    from ip3r.ui.movie_recorder import image_array
    img = QtGui.QImage(3, 2, QtGui.QImage.Format.Format_RGB888)   # 9 bytes a row, padded
    assert img.bytesPerLine() > 9
    for x in range(3):
        for y in range(2):
            img.setPixelColor(x, y, QtGui.QColor(10 * x, 100 + y, 7))
    a = image_array(img)
    assert a.shape == (2, 3, 3)
    assert a[1, 2].tolist() == [20, 101, 7] and a[0, 0].tolist() == [0, 100, 7]
