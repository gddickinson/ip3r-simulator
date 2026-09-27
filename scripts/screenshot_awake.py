"""Keep the smoke test at full speed on macOS.

A GUI app whose window is behind others is throttled by App Nap: the
smoke test was measured at ~30 % of one core, which made every step
several times slower and let a healthy run overrun its hang timer. This
holds a user-initiated, latency-critical activity for the life of the
process (``NSProcessInfo beginActivityWithOptions:reason:``), through the
Objective-C runtime with ctypes, so no PyObjC is needed. Elsewhere it
does nothing.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import sys

__all__ = ["keep_awake"]

# NSActivityUserInitiatedAllowingIdleSystemSleep | NSActivityLatencyCritical
_OPTIONS = 0x00FFFFFF | 0xFF00000000
_held = []


def _send(restype, *argtypes):
    objc = ctypes.cdll.LoadLibrary(ctypes.util.find_library("objc"))
    fn = ctypes.CFUNCTYPE(restype, ctypes.c_void_p, ctypes.c_void_p, *argtypes)
    return objc, fn(ctypes.cast(objc.objc_msgSend, ctypes.c_void_p).value)


def keep_awake(reason: str = "GUI smoke test") -> bool:
    """True when the activity is held (macOS), False elsewhere or on failure."""
    if sys.platform != "darwin":
        return False
    try:
        ctypes.cdll.LoadLibrary(ctypes.util.find_library("Foundation"))
        objc, send = _send(ctypes.c_void_p)
        objc.objc_getClass.restype = ctypes.c_void_p
        objc.sel_registerName.restype = ctypes.c_void_p
        cls, sel = objc.objc_getClass, objc.sel_registerName
        info = send(cls(b"NSProcessInfo"), sel(b"processInfo"))
        _, send_str = _send(ctypes.c_void_p, ctypes.c_char_p)
        text = send_str(cls(b"NSString"), sel(b"stringWithUTF8String:"),
                        reason.encode())
        _, begin = _send(ctypes.c_void_p, ctypes.c_uint64, ctypes.c_void_p)
        token = begin(info, sel(b"beginActivityWithOptions:reason:"), _OPTIONS, text)
        _, retain = _send(ctypes.c_void_p)
        _held.append(retain(token, sel(b"retain")))
        return bool(token)
    except (OSError, AttributeError, TypeError):
        return False

