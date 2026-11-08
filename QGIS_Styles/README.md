# QGIS styles for every raster the workflow writes

One `.qml` per raster, named exactly like the raster (`Zonal_Recovery_Classes.qml`
styles `Zonal_Recovery_Classes.tif`, in Overall and in every barangay folder).

**Apply one style:** in QGIS, right-click the layer → Properties → Symbology →
Style → Load Style… → pick the `.qml` with the same name. Or drag the `.qml` onto
the layer in the Layers panel.

**Apply automatically:** copy a `.qml` next to the `.tif` with the same name; QGIS
uses it whenever that raster is added.

| Style | Raster type | Notes |
|---|---|---|
| `Zonal_Recovery_Classes` (main map), `ITSA_Recovery_Classes`, `Original_*`, `Raw_Pixel_*` | Recovery classes | dark green recovered, light green partial, red not recovered, blue no loss, grey indeterminate |
| `Zonal_RRI_Interpretation_Classes` | Zonal RRI classes | same colours as `Zonal_Recovery_Classes` for the same outcome |
| `RRI_Interpretation_Classes` | Pixel RRI classes | red → orange → yellow → light green → dark green |
| `*_Validity_Status` | RRI validity | green valid, blue no disturbance, grey unresolved |
| `ITSA_Exclusion_Reasons`, `Raw_ITSA_Reasons`, `Zonal_Indeterminate_Reasons` | Why a pixel/zone is Indeterminate | one colour per reason code |
| `Zonal_Disturbance_Evidence` | Disturbance evidence | red event-pair loss, orange modeled loss, purple post-event decline, blue no loss |
| `Zonal_Zone_Size_m` | Zone size used | 70 / 140 / 280 / 560 m |
| `Barangay_Owner` | Barangay of each pixel | labelled with barangay names |
| `ITSA_RRI`, `Zonal_RRI`, `ITSA_Recovery_Fraction` | Index values | −0.5 red … 0 orange … 1 green (back to baseline) … 1.5 dark green |
| `ITSA_Baseline`, `ITSA_Metrics`, `Raw_ITSA_Metrics` | Coherence values (metrics: band 1 = baseline) | viridis 0–0.6 |
| `ITSA_Disturbance_Magnitude`, `ITSA_Final_Deficit`, `ITSA_Recovery_Progress`, `RRI_*_Magnitude` | Signed coherence changes | diverging around 0 |
| `ITSA_Crossing_Days` | Days to return to baseline | viridis 0–120 days |
| `Valid_Data_Coverage`, `ITSA_Gate_Outcomes` | Multi-band yes/no (band 1 shown) | red no, green yes; change band in Symbology |
| `Aligned_Coherence_Observations` | Input coherence stack (band 1 shown) | grey 0–1; change band to view other dates |
| `*_Quality_Flags`, `Spatial_Object_ID`, `Zonal_Zone_ID` | Bit flags / identifiers | greyscale; decode flags with the codebook in the raster metadata |
