"""Defaults for crossfeeding paper figure reproduction"""

MAX_HOURS_DEFAULT = 50

# Simulation / flux conversion constants (PhysiCell_settings.xml)
DFBA_DT_MIN = 0.01
VOXEL_DX_UM = 0.8
VOXEL_VOLUME_UM3 = VOXEL_DX_UM ** 3
CELL_DENSITY_PG_UM3 = 1.04
L_TO_UM3 = 1e15
PG_TO_G = 1e-12
MMOL_TO_FMOL = 1e12
NET_EXPORT_SOLID_FRACTION = 0.75
NET_EXPORT_CELL_DENSITY_G_PER_ML = 1.04
UM3_TO_ML = 1e-12

UM3_PER_L = 1e15
INTERIOR_BOUNDARY_LAYERS = 2

SPECIES_LINE_COLORS = {
    'C. beijerinckii': '#2980B9',
    'M. barkeri': '#D35400',
}

CELL_TYPE_CB = 0
CELL_TYPE_MB = 1

# Bare R_EX_* SBML reaction ids (homogenized with the other use cases); h2/co2/acetate
# are shared between both species, so the same id appears in both lists.
KEY_FLUXES_CB = ['R_EX_glc_D_e', 'R_EX_h2_e', 'R_EX_co2_e', 'R_EX_ac_e']
KEY_FLUXES_MB = ['R_EX_h2_e', 'R_EX_ac_e', 'R_EX_co2_e', 'R_EX_ch4_e']

FLUX_COLORS_CB = ['#f1c232', '#4285f4', '#8e44ad', '#cc0000']
FLUX_COLORS_MB = ['#4285f4', '#cc0000', '#8e44ad', '#2a9d8f']

FLUX_DISPLAY_NAMES = {
    'R_EX_glc_D_e': 'Glucose',
    'R_EX_h2_e': 'H₂',
    'R_EX_co2_e': 'CO₂',
    'R_EX_ac_e': 'Acetate',
    'R_EX_ch4_e': 'CH₄',
}

SUBSTRATE_HEX = {
    'glucose': '#f1c232',
    'acetate': '#cc0000',
}

SUBSTRATE_TS_COLORS = {
    'glucose': SUBSTRATE_HEX['glucose'],
    'methane': '#2a9d8f',
    'H2': '#457b9d',
    'h2': '#457b9d',
    'acetate': SUBSTRATE_HEX['acetate'],
}

SUBSTRATE_DISPLAY_NAMES = {
    'glucose': 'Glucose',
    'methane': 'CH₄',
    'CH4': 'CH₄',
    'ch4': 'CH₄',
    'H2': 'H₂',
    'h2': 'H₂',
    'hydrogen': 'H₂',
    'acetate': 'Acetate',
}
