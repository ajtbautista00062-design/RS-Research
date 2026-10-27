"""Fail-fast spatial settings; physical areas require a projected metric grid."""
import math
import numpy as np
from pyproj import CRS


def validate(c,ref=None):
    if c.OBJECT_HETEROGENEITY_POLICY not in ('quality_only','strict'):
        raise ValueError('Invalid heterogeneity policy.')
    if c.MMU_POLICY not in ('flag_only','defer_low_confidence'):
        raise ValueError('Invalid MMU policy.')
    modes=("pixel_itsa","object_itsa","spatially_regularized_itsa")
    if c.CLASSIFICATION_MODE not in modes:
        raise ValueError(f"CLASSIFICATION_MODE must be one of {modes}")
    if c.CLASSIFICATION_MODE!="pixel_itsa" and not c.ENABLE_SPATIAL_SEGMENTATION:
        raise ValueError("Object modes require ENABLE_SPATIAL_SEGMENTATION.")
    for key in ("MIN_REGION_AREA_M2","TARGET_REGION_AREA_M2","MAX_REGION_AREA_M2",
                "OBJECT_MIN_VALID_PIXELS","OBJECT_MAX_SPATIAL_MAD","MAX_SPATIAL_OBJECTS"):
        value=getattr(c,key)
        if value is None or not np.isfinite(value) or value<=0:
            raise ValueError(f"{key} requires a positive numeric value.")
    if not c.MIN_REGION_AREA_M2<=c.TARGET_REGION_AREA_M2<=c.MAX_REGION_AREA_M2:
        raise ValueError("Require minimum <= target <= maximum object area.")
    if not 0<c.OBJECT_MIN_VALID_FRACTION<=1 or not 0<c.CLASS_CONFIDENCE_THRESHOLD<=1:
        raise ValueError("Invalid object coverage or evidence threshold.")
    if not 0<c.REGULARIZATION_STRENGTH<=1:
        raise ValueError("REGULARIZATION_STRENGTH must be in (0,1].")
    if len(c.PRE_FEATURE_TOLERANCES)!=3 or np.any(np.asarray(c.PRE_FEATURE_TOLERANCES)<=0):
        raise ValueError("Provide three positive PRE feature tolerances.")
    if not isinstance(c.PRE_NEIGHBORHOOD_RADIUS,int) or not 0<=c.PRE_NEIGHBORHOOD_RADIUS<=3:
        raise ValueError("Use a PRE neighborhood radius from 0 to 3.")
    if not c.PRESERVE_HIGH_CONFIDENCE_CLASSES or not c.PRESERVE_NODATA_BOUNDARIES:
        raise ValueError("This implementation requires confidence and NoData protection.")
    if c.VERIFIED_BUILTUP_MASK and not c.BUILTUP_MASK_JUSTIFICATION.strip():
        raise ValueError("Document independent evidence for the built-up mask.")
    if ref is None:
        return None
    crs=CRS.from_user_input(ref.crs)
    if not crs.is_projected or len(crs.axis_info)<2:
        raise ValueError("Object areas require a projected CRS; no degree-based areas.")
    factor=crs.axis_info[0].unit_conversion_factor*crs.axis_info[1].unit_conversion_factor
    area=abs(ref.transform.a*ref.transform.e-ref.transform.b*ref.transform.d)*factor
    if not math.isfinite(area) or area<=0 or area>c.MAX_REGION_AREA_M2:
        raise ValueError("Invalid pixel area or maximum region area below one pixel.")
    if ref.width*ref.height >= np.iinfo('uint32').max:
        raise ValueError("Object identifier grid exceeds uint32 capacity.")
    maximum=int(c.MAX_REGION_AREA_M2//area)
    return dict(pixel_area=area,maximum_pixels=maximum,
                target_pixels=max(1,int(c.TARGET_REGION_AREA_M2//area)),
                minimum_pixels=max(1,math.ceil(c.MIN_REGION_AREA_M2/area)),
                cell_side=max(1,math.isqrt(maximum)))
