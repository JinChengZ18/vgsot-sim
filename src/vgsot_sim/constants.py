from __future__ import annotations

from .configs import PhysicalConstantsConfig

DEFAULT_CONSTANTS = PhysicalConstantsConfig()
_active_constants = PhysicalConstantsConfig()


def get_constants() -> PhysicalConstantsConfig:
    return _active_constants


def set_constants(constants: PhysicalConstantsConfig) -> PhysicalConstantsConfig:
    global _active_constants
    _active_constants = constants
    return _active_constants


def reset_constants() -> PhysicalConstantsConfig:
    return set_constants(PhysicalConstantsConfig())


def __getattr__(name: str):
    constants = get_constants()
    if hasattr(constants, name):
        return getattr(constants, name)
    raise AttributeError(f"module 'constants' has no attribute {name!r}")


__all__ = [
    'PhysicalConstantsConfig',
    'DEFAULT_CONSTANTS',
    'get_constants',
    'set_constants',
    'reset_constants',
]
