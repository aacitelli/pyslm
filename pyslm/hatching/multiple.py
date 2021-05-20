from pyslm.hatching.hatching import Hatcher
from typing import Any 

# TODO: Probably cleaner to take a list of 2-tuples in. Not a big deal though.
# TODO: Figure out actual types here instead of being lazy and doing Any. Just for type hints and doesn't affect runtime, but this is a library and supposed to be readable and maintainable.
# TODO: Finish documentation 
def hatch_multiple(default_hatcher: Hatcher, hatcher_area_pairs: tuple[Hatcher, list], geometry: Any) -> None: 
    """Given an array of pairs of hatchers and associated areas for them to hatch, returns 

    :param geometry: 
    :type geometry: The polygon to hatch, provided as whatever the return value of Part.getVectorSlice() is. 
        Note that we don't really care what it is, as we simply propagate it to hatch() calls.
    :param hatchers: List of 2-Long Lists, where the first element of each is an initialized hatcher (i.e. all parameters set) 
        and the second is a four-long array specifying a [min-x, min-y, max-x, max-y] to trim the resulting hatches to.
        Note that a default hatcher will be initialized, run, and trimmed such that it includes everything except for the specified areas.
    :type hatchers: list
    """

    # 1. Iterate through Hatchers, getting hatching result then trimming result to associated area 
    for hatcher_area_pair in hatcher_area_pairs:
        layer = hatcher_area_pair[0].hatch(geometry)
        trim_layer_to_inside_area(layer, hatcher_area_pair[1])

    # 2. Run default hatcher and trim out any areas that 
    default_hatches = default_hatcher.hatch()
    trim_layer_to_outside_areas(layer, hatcher_area_pairs)

    # 3. Concatenate all vectors into one representation and return 

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