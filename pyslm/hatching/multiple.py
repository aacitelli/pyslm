# TODO: Organize imports

from numpy import isin
from pyslm.hatching.hatching import Hatcher
from pyslm.geometry.geometry import Layer, HatchGeometry, ContourGeometry
from typing import Any, Union

# Used for line intersection
import shapely
from shapely.geometry import LineString, Point, point
import numpy as np
import numpy.typing as npt

# TODO: Remove
from typeguard import typechecked

# TODO: Function does a lot of number crunching; convert everything to numpy operations, if they aren't yet. Cleanest way to do this is likely to run a profiler.

# TODO: Probably cleaner to take a list of 2-tuples in. Not a big deal though.
# TODO: Figure out actual types here instead of being lazy and doing Any. Just for type hints and doesn't affect runtime, but this is a library and supposed to be readable and maintainable.


@typechecked  # TODO: Remove
def hatch_multiple(hatchers: list[Hatcher], areas: npt.ArrayLike, default_hatcher: Hatcher, boundary: npt.ArrayLike, z: float) -> Layer:
    """Runs hatchers on their associated areas and the default hatcher on anything not in an area, and returns the
    concatenation of all the hatches. Ordered by the order they are passed-in. 

    :param hatchers: A `list` of Hatcher objects initialized with all their parameters.
    :type hatchers: list
    :param areas: A `list` of areas, each a `list` in [min-x, min-y, max-x, max-y] form 
    :type areas: list
    :param default_hatcher: A Hatcher object initialized with all its parameters that should be used at any point not in an area.
    :type default_hatcher: Hatcher
    :param boundary: The layer slice returned by `Part.getVectorSlice(z-coordinate)`; generally an ndarray of float points representing consecutive points
    :type boundary: Any
    :param height: The current z-axis height of the given layer (needed to trim areas); in mm 
    :type height: float 
    """

    # 1. Iterate through Hatchers, getting hatching result then trimming result to associated area
    trimmed_layers = []
    for i in range(len(hatchers)):

        # If this hatcher shouldn't apply to this z-coordinate, skip it
        min_z, max_z = areas[i][2], areas[i][5]
        if z > max_z or z < min_z:
            continue

        # Run hatcher on contour boundary
        layer = hatchers[i].hatch(boundary)

        # Only keep contour of the first Hatcher (to avoid duplicates)
        if i != 0:
            layer.geometry = [
                x for x in layer.geometry if isinstance(x, HatchGeometry)]

        # Trims any instances of HatchGeometry to the given area
        trimmed_layers.append(trim_layer_to_inside_area(layer, areas[i]))

    # Edge Case: None of the areas applied at this height, so we should hatch with the default and return it
    if not len(trimmed_layers):
        return default_hatcher.hatch(boundary)

    # 2. Concatenate all trimmed vectors into one `Layer` object
    end_layer = trimmed_layers[0]
    trimmed_layers = trimmed_layers[1:]
    for layer in trimmed_layers:
        for geometry in layer.geometry:
            end_layer.geometry.append(geometry)

    # 3. TODO: Run default hatcher, trim out applicable areas, then add the geometries
    # At least one custom hatcher was used; TODO: Hatch with default then cut out any applicable areas
    # default_hatches = default_hatcher.hatch(boundary)
    # trim_layer_to_outside_areas(layer, areas)

    return end_layer

# TODO: Function works with a high amount of data, make sure these are all NumPy operations
# TODO: Figure out type hints


@typechecked  # TODO: Remove
def trim_layer_to_inside_area(layer: Layer, area: npt.ArrayLike) -> Layer:
    """Trims any HatchGeometry instances on the layer to the provided area 

    :param layer: The layer object, generally with HatchGeometry instances in it 
    :type layer: pyslm.
    :param area: Area to trim the specified hatches to; organized 
    :type area: Any
    :return: The provided Layer, with hatches trimmed to the given area. Returned in the same order.
    :rtype: Any
    """

    # Trim any HatchGeometry instances on the layer to the provided area
    # NOTE: If this function is called, we know the given Hatcher was valid for this z-coordinate
    end_geometries = []
    for geometry in layer.geometry:
        if isinstance(geometry, ContourGeometry):  # Don't modify contour geometry
            end_geometries.append(geometry)
            continue

        # We don't know how big the trimmed version will be, so more efficient to use Python arr then convert to np at end; see https://stackoverflow.com/a/10122262/6402548
        trimmed_hatches = []
        for i in range(0, len(geometry.coords), 2):
            segment = LineString([geometry.coords[i], geometry.coords[i + 1]])
            trimmed_segment = trim_segment_to_box(segment, area)
            if trimmed_segment != None:
                trimmed_hatches.append(
                    [trimmed_segment.xy[0][0], trimmed_segment.xy[1][0]])
                # xy is represented as [x1, x2],[y1, y2], so we essentially do a hstack
                trimmed_hatches.append(
                    [trimmed_segment.xy[0][1], trimmed_segment.xy[1][1]])
        if len(trimmed_hatches):
            geometry.coords = np.array(trimmed_hatches)
            end_geometries.append(geometry)
    layer.geometry = end_geometries

    return layer

# TODO: Function works with a high amount of data, make sure these are all NumPy operations


@typechecked  # TODO: Remove
def trim_layer_to_outside_areas(layer: Layer, area: npt.ArrayLike) -> Layer:
    """[summary]

    :param hatches: List of hatches to trim to OUTSIDE OF the specified area.
    :type hatches: Any
    :param area: Full list of <Hatcher, area specification> pairs to trim the specified hatches to OUTSIDE OF.
    :type area: Any
    :return: The provided layer with (1) normal contours and (2) hatches trimmed to OUTSIDE OF the given area. Returned in the same order.
    :rtype: Layer
    """

    # Other stuff still isn't working, so this is on the backburner
    raise NotImplementedError()

    # Trim any HatchGeometry instances on the layer to the provided area
    # NOTE: If this function is called, we know the given Hatcher was valid for this z-coordinate
    """
    for geometry in layer.geometry:
        if isinstance(geometry, ContourGeometry):
            continue
        for i in range(0, len(geometry.coords), 2):
            segment = LineString([geometry.coords[i], geometry.coords[i + 1]])
            trimmed_segment = trim_segment_to_box(segment, area)
            if trimmed_segment != None:
                geometry.coords[i], geometry.coords[i + 1] = trimmed_segment
    """

# Largely pulled from a separate implementation I previously did; Python/Shapely makes some stuff more concise, though
# See https://github.com/aacitelli/OASIS_IslandScanning/blob/fce5e3fb72c71dc22beb887747639b3969390c6b/Source_Code_Ohio_State_CDME/genScan/ScanPath.cpp#L924
@typechecked  # TODO: Remove
def trim_segment_to_box(segment: LineString, area: np.ndarray) -> Union[LineString, None]:
    """Trims a given Shapely `LineString` object to the given area, and returns 

    :param segment: The `LineString` object to trim to the given area.
    :type segment: `shapely.geometry.LineString`
    :param area: The area to trim the given `LineString` object to.
    :type area: list
    :return: Either:
    - A Shapely `LineString` object, f a subsegment of the line is inside the provided area
    - `None`, if the given segment has zero overlap with the provided area.
    :rtype: Union[`shapely.geometry.LineString`, Any]
    """

    # Misc. Helper Functions
    @typechecked  # TODO: Remove
    def point_in_bbox(v: Point, bl: Point, tr: Point):
        return v.x >= bl.x and v.x <= tr.x and v.y >= bl.y and v.y <= tr.y

    @typechecked  # TODO: Remove
    def vertices_equal(v1: Point, v2: Point):
        return abs(v1.x - v2.x) <= .00000001 and abs(v1.y - v2.y) <= .00000001

    # Renaming and putting in formats more readable/usable
    v1, v2 = Point(segment.xy[0][0], segment.xy[1][0]), Point(segment.xy[0][1], segment.xy[1][1])
    min_x, min_y, max_x, max_y = area[0], area[1], area[3], area[4]
    top_left, top_right, bottom_left, bottom_right = Point(min_x, max_y), Point(max_x, max_y), Point(min_x, min_y), Point(max_x, min_y)
    top_segment, right_segment = LineString([top_left, top_right]), LineString([top_right, bottom_right])
    bottom_segment, left_segment = LineString([bottom_right, bottom_left]), LineString([bottom_left, top_left])

    # Get intersection points on each side
    top_intersection = segment.intersection(top_segment)
    right_intersection = segment.intersection(right_segment)
    bottom_intersection = segment.intersection(bottom_segment)
    left_intersection = segment.intersection(left_segment)

    # LineString.intersection(LineString) returns the following general cases (documenting b/c the documentation sucks)
    # If intersection, returned as a single shapely.geometry.point.Point object
    # If no intersection or an intersection would only occur if segment(s) were extended, len(intersection.coords) will be zero
    # If they are colinear, LineString with xy being an array of the points

    # Non-empty LineString return => Segment was colinear with one of the edges
    # By personal convention, only return Top/Left borders to avoid duplicate vectors on area boundaries
    if isinstance(top_intersection, LineString) and len(top_intersection.coords) != 0:
        return LineString([top_intersection.coords[0], top_intersection.coords[1]])
    if isinstance(left_intersection, LineString) and len(left_intersection.coords) != 0:
        return LineString([left_intersection.coords[0], left_intersection.coords[1]])
    if isinstance(right_intersection, LineString) and len(right_intersection.coords) != 0 or \
            isinstance(bottom_intersection, LineString) and len(bottom_intersection.coords) != 0:
        return None

    # Trim sides where there was no intersection (which will be an empty LineString)
    intersections = [top_intersection, right_intersection,
                     bottom_intersection, left_intersection]
    intersections = [x for x in intersections if isinstance(x, Point)]

    # Trim duplicate intersections (only occurs when segment ends at a corner)
    intersections_without_duplicates = []
    for intersection in intersections:
        for unique_intersection in intersections_without_duplicates:
            if vertices_equal(intersection, unique_intersection):  # Current coordinates already added
                break
        intersections_without_duplicates.append(
            intersection)  # Gets through loop => It's unique
    intersections = intersections_without_duplicates

    # Zero distinct intersections => Segment is either completely inside or outside, only necessary to test one point
    if len(intersections) == 0:
        if v1.x < top_left.x or v1.x > top_right.x or \
                v1.y < bottom_left.y or v1.y > top_left.y:
            return None
        # Could probably just return `segment`, but don't want to risk weird reference behavior
        return LineString([v1, v2])

    # One distinct intersection => Subsegment is the intersection point and whichever endpoint is inside (or on) the boundary
    elif len(intersections) == 1:
        if point_in_bbox(v1, bottom_left, top_right) and not vertices_equal(v1, intersections[0]):
            return LineString([v1, intersections[0]])
        elif point_in_bbox(v2, bottom_left, top_right) and not vertices_equal(v2, intersections[0]):
            return LineString([v2, intersections[0]])
        return None  # Segments that end exactly on boundary... won't bother adding a zero-length vector 

    # Two distinct intersections => Subsegment
    elif len(intersections) == 2:
        return LineString([intersections[0], intersections[1]])

    else:
        raise Exception(
            "Got invalid number of distinct intersections: {}".format(len(intersections)))
