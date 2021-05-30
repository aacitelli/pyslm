from pyslm.hatching.hatching import Hatcher
from typing import Any 

# TODO: Probably cleaner to take a list of 2-tuples in. Not a big deal though.
# TODO: Figure out actual types here instead of being lazy and doing Any. Just for type hints and doesn't affect runtime, but this is a library and supposed to be readable and maintainable.
# TODO: Finish documentation 
def hatch_multiple(hatchers: list, areas: list, default_hatcher: Hatcher, boundary: Any) -> None: 
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
    """

    # 1. Iterate through Hatchers, getting hatching result then trimming result to associated area 
    hatches = []
    for i in range(len(hatchers)):
        layer = hatchers[i].hatch(boundary)
        trimmed_layer = trim_layer_to_inside_area(layer, areas[i])
        hatches.append(trimmed_layer)

    # 2. Run default hatcher and trim out any areas that 
    default_hatches = default_hatcher.hatch(boundary)
    trim_layer_to_outside_areas(layer, areas)

    # 3. TODO: Concatenate all vectors into one representation and return 


# TODO: Function works with a high amount of data, make sure these are all NumPy operations 
# TODO: Figure out type hints
def trim_layer_to_inside_area(layer: Any, area: Any) -> Any:
    """[summary]

    :param hatches: List of hatches to trim to the specified area.
    :type hatches: Any
    :param area: Area to trim the specified hatches to.
    :type area: Any
    :return: The provided list of hatches trimmed to the given area. Returned in the same order.
    :rtype: Any
    """

    raise NotImplementedError("TODO: Implement!")

# TODO: Function works with a high amount of data, make sure these are all NumPy operations 
# TODO: Figure out type hints
def trim_layer_to_outside_areas(layer: Any, area: Any) -> Any:
    """[summary]

    :param hatches: List of hatches to trim to OUTSIDE OF the specified area.
    :type hatches: Any
    :param area: Full list of <Hatcher, area specification> pairs to trim the specified hatches to OUTSIDE OF.
    :type area: Any
    :return: The provided list of hatches trimmed to OUTSIDE OF the given area. Returned in the same order.
    :rtype: Any
    """

    raise NotImplementedError("TODO: Implement!")