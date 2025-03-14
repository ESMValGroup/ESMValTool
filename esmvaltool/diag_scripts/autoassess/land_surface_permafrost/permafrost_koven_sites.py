"""
Module containing coordinates of observed sites from Koven paper.

Written in format required by Iris interpolator:
    (('latitude',  <lats_tuple>),('longitude', <lons_tuple>))
"""

import numpy as np

lats = np.array(
    [
        71.320575,
        64.9076,
        64.86781,
        68.0691666666667,
        70.161283,
        70.161283,
        70.1612833333333,
        69.6551333333333,
        69.6741402,
        69.6741402,
        69.65513,
        69.65513,
        68.4776666666667,
        70.31516,
        70.31516,
        69.14664082,
        69.14664082,
        69.239,
        68.6983,
        67.0134,
        65.3158,
        68.2875,
        68.2903,
        65.3134,
        65.3146,
        56.761,
        69.4283,
        69.43303,
        66.9357,
        66.9381,
        64.86936667,
        64.86611667,
        64.86751667,
        64.86936,
        64.867516,
        64.867516,
        64.8669,
        70.3744666666667,
    ]
)

lons = np.array(
    [
        -156.6493305,
        -163.674483,
        -147.78486111,
        -149.580333333333,
        -148.4653,
        -148.4653,
        -148.4653,
        -148.722016666667,
        -148.72076632,
        -148.72076632,
        -148.722016,
        -148.722016,
        -149.501666666667,
        -147.99316,
        -147.99316,
        -148.8483008,
        -148.8483008,
        -51.0623,
        -149.35181,
        -50.7175,
        72.8745,
        54.5026,
        54.5026,
        72.8872,
        72.8629,
        118.1901,
        -148.70015,
        -148.67385,
        -53.6416,
        -53.5957,
        -147.8608,
        -147.8567667,
        -147.8588333,
        -147.8608,
        -147.85676,
        -147.85883,
        -147.858383,
        -148.552166666667,
    ]
)

site_points = (("latitude", lats), ("longitude", lons))
