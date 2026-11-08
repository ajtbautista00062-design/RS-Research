<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <colorPalette>
        <paletteEntry value="0" color="#000000" alpha="0" label="0 Outside study area / structurally unavailable"/>
        <paletteEntry value="1" color="#1a9641" alpha="255" label="1 Valid RRI"/>
        <paletteEntry value="2" color="#2c7bb6" alpha="255" label="2 No detectable coherence disturbance"/>
        <paletteEntry value="3" color="#fdae61" alpha="255" label="3 Insufficient observations"/>
        <paletteEntry value="4" color="#f46d43" alpha="255" label="4 Invalid or near-zero disturbance denominator"/>
        <paletteEntry value="5" color="#a50026" alpha="255" label="5 Invalid fitted model result"/>
        <paletteEntry value="6" color="#bdbdbd" alpha="255" label="6 Unresolved or contradictory temporal evidence"/>
        <paletteEntry value="255" color="#000000" alpha="0" label="255 NoData"/>
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" invertColors="0" colorizeOn="0" colorizeRed="255" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
