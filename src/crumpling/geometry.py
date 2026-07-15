"""Numerical geometry helpers."""

import math

import numpy as np


def maximum_pairwise_distance(points: np.ndarray) -> float:
    """Return the largest Euclidean separation without an O(n²) matrix."""

    if len(points) < 2:
        return 0.0
    largest_squared = 0.0
    for index in range(len(points) - 1):
        differences = points[index + 1 :] - points[index]
        squared = float(np.max(np.einsum("ij,ij->i", differences, differences)))
        largest_squared = max(largest_squared, squared)
    return math.sqrt(largest_squared)


def radius_of_gyration(points: np.ndarray) -> float:
    """Return root-mean-square distance from the centroid."""

    if len(points) == 0:
        raise ValueError("radius of gyration requires at least one point")
    differences = points - points.mean(axis=0)
    return math.sqrt(float(np.einsum("ij,ij->", differences, differences)) / len(points))


def planar_inertia_eigenvalues(points: np.ndarray) -> tuple[float, float]:
    """Return the principal moments of the x-y point cloud."""

    if len(points) <= 1:
        return 0.0, 0.0
    center = points.mean(axis=0)
    dx = center[0] - points[:, 0]
    dy = center[1] - points[:, 1]
    i_xx = float(np.sum(dy * dy))
    i_yy = float(np.sum(dx * dx))
    i_xy = float(-np.sum(dx * dy))
    discriminant = math.sqrt((i_xx - i_yy) ** 2 + 4.0 * i_xy**2)
    return (
        0.5 * (i_xx + i_yy + discriminant),
        0.5 * (i_xx + i_yy - discriminant),
    )
