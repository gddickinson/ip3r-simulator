"""File → Record movie…'s defaults (the dialog starts here; every one can
be changed there). Drawing only: no number here enters a calculation.
Imported into ``parameter_table.P``.
"""

from param_entry import entry as _p

MOVIE = [
    _p("movie.fps", "Movie frame rate", 20, "frames/s", "convention", "display",
       "convention", "Playback rate written into a recorded GIF or MP4.",
       "A GIF holds delays in hundredths of a second, and 20 frames/s (50 ms) "
       "is the nearest smooth rate it holds exactly", 1, 60),
    _p("movie.turntable_frames", "Frames per camera turn", 72, "frames",
       "convention", "display", "convention", "A turntable movie turns the "
       "camera 360 degrees in this many frames (5 degrees a frame).",
       "At 20 frames/s one turn takes 3.6 s", 2, 720),
    _p("movie.mode_frames", "Frames per mode cycle", 48, "frames", "convention",
       "display", "convention", "A normal-mode movie shows one sinusoidal "
       "cycle of the animating mode in this many frames.",
       "At 20 frames/s one cycle takes 2.4 s; the mode's time scale is "
       "illustrative (an elastic network gives shapes, not rates)", 2, 720),
    _p("movie.hold_frames", "Frames held at each end state", 10, "frames",
       "convention", "display", "convention", "A transition movie holds the "
       "start and the end deposit this many frames each before turning back.",
       "Half a second at 20 frames/s, long enough to see each state", 0, 120),
    _p("movie.webp_quality", "WebP quality", 75, "", "convention", "display",
       "convention", "Lossy quality (0-100) of a recorded animated WebP.",
       "At 75 an 8TKG -> 8TKF morph is 1.6 MB against 4.1 MB as a GIF, with no "
       "visible loss at the recorded size (2026-09-28)", 1, 100),
    _p("movie.max_width", "Movie width", 720, "px", "convention", "display",
       "convention", "Frames are scaled down (never up) to at most this width.",
       "A 720 px GIF of 60-80 frames stays within a few MB for a README", 64, 3840),
]
