<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <colorPalette>
        <paletteEntry value="0" color="#000000" alpha="0" label="0 NoData / RRI not valid"/>
        <paletteEntry value="1" color="#d7191c" alpha="255" label="1 Further coherence degradation"/>
        <paletteEntry value="2" color="#d7191c" alpha="255" label="2 No measurable recovery"/>
        <paletteEntry value="3" color="#a6d96a" alpha="255" label="3 Partial coherence recovery"/>
        <paletteEntry value="4" color="#1a9641" alpha="255" label="4 Return to modeled baseline"/>
        <paletteEntry value="5" color="#1a9641" alpha="255" label="5 Coherence exceeding modeled baseline"/>
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" invertColors="0" colorizeOn="0" colorizeRed="255" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
