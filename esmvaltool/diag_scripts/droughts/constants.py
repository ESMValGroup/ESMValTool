"""Constants that might be useful by any drought diagnostic."""

from iris.coords import AuxCoord

# fmt: off
DENSITY = AuxCoord(1000, long_name="density", units="kg m-3")

FNAME_FORMAT = "{project}_{reference_dataset}_{mip}_{exp}_{ensemble}_"
"{short_name}_{start_year}-{end_year}"
CMIP6_FNAME = "{project}_{dataset}_{mip}_{exp}_{ensemble}_{short_name}_"
"{grid}_{start_year}-{end_year}"
OBS_FNAME = "{project}_{dataset}_{type}_{version}_{mip}_{short_name}_"
"{start_year}-{end_year}"

CONTINENTAL_REGIONS = {
    "Global": ["GLO"],
    "North America": ["GIC", "NWN", "NEN", "WNA", "CNA", "ENA"],
    "Central America": ["NCA", "SCA", "CAR"],
    "Southern America": ["NWS", "NSA", "NES", "SAM", "SWS", "SES", "SSA"],
    "Europe": ["NEU", "WCE", "EEU", "MED"],
    "Africa": ["SAH", "WAF", "CAF", "NEAF", "SEAF", "WSAF", "ESAF", "MDG"],
    "Asia": ["RAR", "WSB", "ESB", "RFE", "WCA", "ECA", "TIB", "EAS", "ARP",
             "SAS", "SEA"],
    "Australia": ["NAU", "CAU", "EAU", "SAU", "NZ", "WAN", "EAN"],
}

HEX_POSITIONS = {
    "NWN": [2, 0], "NEN": [4, 0], "GIC": [6.5, -0.5], "NEU": [14, 0],
    "RAR": [20, 0], "WNA": [1, 1], "CNA": [3, 1], "ENA": [5, 1],
    "WCE": [13, 1], "EEU": [15, 1], "WSB": [17, 1], "ESB": [19, 1],
    "RFE": [21, 1], "NCA": [2, 2], "MED": [14, 2], "WCA": [16, 2],
    "ECA": [18, 2], "TIB": [20, 2], "EAS": [22, 2], "SCA": [3, 3],
    "SAH": [13, 3], "ARP": [15, 3], "SAS": [19, 3], "SEA": [23, 3],
    "NWS": [6, 4], "NSA": [8, 4], "WAF": [12, 4], "CAF": [14, 4],
    "NEAF": [16, 4], "NAU": [24.5, 4.3], "SAM": [7, 5], "NES": [9, 5],
    "WSAF": [13, 5], "SEAF": [15, 5], "MDG": [17.5, 5.3],
    "CAU": [23.5, 5.3], "EAU": [25.5, 5.3], "SWS": [6, 6], "SES": [8, 6],
    "ESAF": [14, 6], "SAU": [24.5, 6.3], "NZ": [27, 6.5], "SSA": [7, 7],
}

MAINLY_SEA_HEX_POSITIONS = {
    "PAC": [27.5, 3.3],  # Pacific
    "CAR": [5, 3],  # Caribbean
}

INDEX_META = {
    "CDD": {
        "long_name": "Conscutive Dry Days",
        "short_name": "CDD",
        "units": "days",
        "standard_name": "consecutive_dry_days",
    },
    "PDSI": {
        "long_name": "Palmer Drought Severity Index",
        "showrt_name": "PDSI",
        "units": "1",
        "standard_name": "palmer_drought_severity_index",
    },
    "SCPDSI": {
        "long_name": "Self-calibrated Palmer Drought Severity Index",
        "short_name": "scPDSI",
        "units": "1",
        "standard_name": "self_calibrated_palmer_drought_severity_index",
    },
    "SPI": {
        "long_name": "Standardized Precipitation Index",
        "short_name": "SPI",
        "units": "1",
        "standard_name": "standardized_precipitation_index",
    },
    "SPEI": {
        "long_name": "Standardized Precipitation Evapotranspiration Index",
        "short_name": "SPEI",
        "units": "1",
        "standard_name": "standardized_precipitation_evapotranspiration_index",
    },
}
# fmt: on
