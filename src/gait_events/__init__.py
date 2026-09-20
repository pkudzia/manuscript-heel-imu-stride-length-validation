"""Gait event detection.

The pipeline detects heel strikes and toe-offs with an adaptive-threshold
detector adapted from Salarian et al. (2004), with signal-derived
thresholds in the spirit of Greene et al. (2010). Adaptive thresholds
matter here because gait speed, and therefore gyroscope amplitude, varies
widely across older-adult participants and cognitive-load conditions.
"""

from .adaptive import AdaptiveThresholdDetector

__all__ = ["AdaptiveThresholdDetector"]
