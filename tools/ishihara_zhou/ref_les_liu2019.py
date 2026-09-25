"""LES reference for the 2D ridge: Liu, Diao & Ishihara, Renewable Energy 136
(2019) 968-992, digitized from the 150 dpi author PDF.

Smooth ridge (2Ds): smooth wall. Rough ridge (2Dr): 5 mm canopy drag layer
(artificial grass). Their experiment is the Ishihara-Hibi measurement of the
same ridge, not the 2025 Data in Brief files.

U/Uref is read from Fig. 8c at 5 mm (smooth) or 7 mm (rough) above the local
surface, the heights of the lowest 2025 tunnel points. Uref is their flat
terrain U at 4 h; it is converted with the 2025 approach-flow U at 160 mm
(mmc1 / mmc2). Digitizing error is about 0.05-0.1 Uref (0.3-0.5 m/s).

The reversed-flow extent at the same height is from the bubble outlines of
Fig. 11 (top of the U < 0 region), x/h from the crest.
"""

UREF = {"smooth": 5.3923, "rough": 5.2973}
STATIONS_XH = (-2.5, -1.25, 0.0, 1.25, 2.5, 3.75, 5.0, 6.25)
U_OVER_UREF = {
    "smooth": (0.40, 0.63, 1.05, -0.05, -0.15, -0.08, 0.03, 0.18),
    "rough": (0.15, 0.38, 0.69, -0.05, -0.02, -0.11, -0.06, 0.02),
}
REVERSED_XH = {"smooth": (1.12, 4.50), "rough": (0.93, 5.97)}
H = 40.0  # ridge height in m at full scale
