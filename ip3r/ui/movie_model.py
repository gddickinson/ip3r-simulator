"""Movies of the 3-D view, Qt-free: what each frame shows, and the file
the frames become.

A movie is a *plan* of frames, each a camera angle and a scene state, so the
same plan gives the same movie at any machine speed: nothing is sampled
from the live animation's clock, which drops frames when the machine is
busy. Three things can move:

- ``turntable``: the camera turns once about the screen's vertical axis
  (the four-fold axis in the side view), the structure still;
- ``mode``: the animated normal mode through one cycle,
  ``xyz = base + sin(phase) · displacement`` (``modes_panel.mode_frame``);
- ``transition``: the morph's solved frames there and back, held at each
  end state. Only solved frames are shown, never a blend of two (a blend of
  two restrained frames is not itself restrained).

A mode or a transition can also turn the camera once while it plays.

The frames are loops: a turntable or a mode cycle stops one step short of
where it began, so the last frame leads straight back into the first.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["SOURCES", "SOURCE_LABELS", "FORMATS", "MovieSpec", "Frame", "plan",
           "pingpong", "fit_size", "formats_available", "gif_delay_ms", "delay_ms",
           "write_gif", "write_webp", "write_mp4",
           "write_movie"]

SOURCES = ("turntable", "mode", "transition")
SOURCE_LABELS = {"turntable": "Turntable (the camera turns once)",
                 "mode": "Normal mode (the one animating)",
                 "transition": "Transition (the built morph, there and back)"}
#: Format → file suffix. GIF plays everywhere; animated WebP keeps full
#: colour at about a third of the size (browsers and GitHub play it); MP4 is
#: for slides and players.
FORMATS = {"gif": ".gif", "webp": ".webp", "mp4": ".mp4"}


@dataclass(frozen=True)
class MovieSpec:
    source: str
    #: Frames per turn (turntable) or per cycle (mode); a transition shows
    #: its own solved frames and ignores this.
    frames: int
    fps: float
    #: Also turn the camera once over the movie (always, for a turntable).
    turn: bool = False
    #: Frames each end state of a transition is held.
    hold: int = 0
    hud: bool = True
    max_width: int = 720
    fmt: str = "gif"

    def __post_init__(self) -> None:
        if self.source not in SOURCES:
            raise ValueError(f"unknown movie source {self.source!r}; known: {SOURCES}")
        if self.fmt not in FORMATS:
            raise ValueError(f"unknown movie format {self.fmt!r}; known: {tuple(FORMATS)}")
        if self.frames < 2 or self.fps <= 0 or self.max_width < 16 or self.hold < 0:
            raise ValueError("a movie needs ≥ 2 frames, a positive rate, a width "
                             "of at least 16 px and a non-negative hold")


@dataclass(frozen=True)
class Frame:
    #: Camera turn from the starting view, degrees.
    angle: float
    #: Mode phase (radians), transition frame index, or None (still).
    state: float | int | None


def pingpong(n_solved: int, hold: int = 0) -> list[int]:
    """Solved frames 0 → n−1 → 1, each end held ``hold`` extra frames; the
    last index leads back into 0, so the movie loops without a repeat."""
    if n_solved < 2:
        raise ValueError("a transition needs at least two solved frames")
    up = list(range(n_solved))
    down = list(range(n_solved - 2, 0, -1))
    return [0] * hold + up + [n_solved - 1] * hold + down


def plan(spec: MovieSpec, n_solved: int | None = None) -> list[Frame]:
    """Every frame of the movie ``spec`` describes. A transition needs the
    number of solved frames of the morph on screen."""
    if spec.source == "turntable":
        states: list = [None] * spec.frames
    elif spec.source == "mode":
        states = list(2.0 * np.pi * np.arange(spec.frames) / spec.frames)
    else:
        if n_solved is None:
            raise ValueError("a transition movie needs the morph's frame count")
        states = pingpong(n_solved, spec.hold)
    n = len(states)
    turn = spec.turn or spec.source == "turntable"
    return [Frame(360.0 * k / n if turn else 0.0, s) for k, s in enumerate(states)]


def fit_size(width: int, height: int, max_width: int, even: bool = False) -> tuple[int, int]:
    """(w, h) no wider than ``max_width``, aspect kept; ``even`` rounds both
    down to even numbers (H.264 in yuv420p needs them)."""
    scale = min(1.0, max_width / width)
    w, h = max(1, int(round(width * scale))), max(1, int(round(height * scale)))
    if even:
        w, h = max(2, w - w % 2), max(2, h - h % 2)
    return w, h


def formats_available() -> dict[str, str | None]:
    """Format → None when it can be written, else why not."""
    from PIL import features
    out: dict[str, str | None] = {"gif": None, "webp": None if features.check("webp")
                                  else "this Pillow was built without WebP"}
    try:
        import imageio_ffmpeg  # noqa: F401
        out["mp4"] = None
    except ImportError:
        out["mp4"] = "MP4 needs the imageio-ffmpeg package (pip install imageio-ffmpeg)"
    return out


def _check(frames: list[np.ndarray]) -> tuple[int, int]:
    if not frames:
        raise ValueError("no frames to write")
    h, w = frames[0].shape[:2]
    for f in frames:
        if f.shape != (h, w, 3) or f.dtype != np.uint8:
            raise ValueError("every frame must be the same H x W x 3 uint8 image")
    return w, h


def gif_delay_ms(fps: float) -> int:
    """GIF stores a frame's delay in hundredths of a second, so the rate is
    rounded to the nearest one it can hold (20 frames/s is exact; 15 would
    play at 16.7)."""
    return 10 * max(1, int(round(100.0 / fps)))


def delay_ms(fps: float, fmt: str) -> int:
    """Each frame's delay as the format stores it (GIF in centiseconds)."""
    return gif_delay_ms(fps) if fmt == "gif" else int(round(1000.0 / fps))


def write_gif(path, frames: list[np.ndarray], fps: float) -> Path:
    """An endlessly looping GIF. One palette is built from frames across the
    whole movie and used for all of them, so a colour does not flicker
    between frames as it would with a palette per frame. Identical
    consecutive frames (a transition's held end states) are stored once with
    their delays added, so the file is smaller and the timing unchanged."""
    from PIL import Image
    _check(frames)
    step = max(1, len(frames) // 8)
    sample = np.concatenate(frames[::step], axis=0)
    palette = Image.fromarray(sample).quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    images = [Image.fromarray(f).quantize(palette=palette, dither=Image.Dither.NONE)
              for f in frames]
    path = Path(path)
    images[0].save(path, save_all=True, append_images=images[1:],
                   duration=gif_delay_ms(fps), loop=0, disposal=1)
    return path


def write_webp(path, frames: list[np.ndarray], fps: float, quality: int = 75) -> Path:
    """An endlessly looping animated WebP in full colour (lossy at
    ``quality``, the ``movie.webp_quality`` parameter by default in the app)."""
    from PIL import Image
    _check(frames)
    reason = formats_available()["webp"]
    if reason:
        raise RuntimeError(reason)
    images = [Image.fromarray(f) for f in frames]
    path = Path(path)
    images[0].save(path, save_all=True, append_images=images[1:],
                   duration=delay_ms(fps, "webp"), loop=0, quality=int(quality), method=6)
    return path


def write_mp4(path, frames: list[np.ndarray], fps: float) -> Path:
    """H.264 in yuv420p (plays in browsers and on GitHub); odd sizes are
    cropped by one pixel, as the codec requires even ones."""
    reason = formats_available()["mp4"]
    if reason:
        raise RuntimeError(reason)
    import imageio_ffmpeg
    w, h = _check(frames)
    w, h = w - w % 2, h - h % 2
    path = Path(path)
    writer = imageio_ffmpeg.write_frames(str(path), (w, h), fps=fps, codec="libx264",
                                         pix_fmt_out="yuv420p", quality=8,
                                         macro_block_size=1)
    writer.send(None)
    try:
        for f in frames:
            writer.send(np.ascontiguousarray(f[:h, :w]).tobytes())
    finally:
        writer.close()
    return path


def write_movie(path, frames: list[np.ndarray], spec: MovieSpec) -> Path:
    """The frames in ``spec``'s format (the suffix is added if missing)."""
    path = Path(path)
    if path.suffix.lower() != FORMATS[spec.fmt]:
        path = path.with_suffix(FORMATS[spec.fmt])
    if spec.fmt == "webp":
        from ..parameters import PARAMETERS as _P
        return write_webp(path, frames, spec.fps, int(_P.value("movie.webp_quality")))
    writer = write_gif if spec.fmt == "gif" else write_mp4
    return writer(path, frames, spec.fps)
