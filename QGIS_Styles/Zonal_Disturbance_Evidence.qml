<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <colorPalette>
        <paletteEntry value="0" color="#bdbdbd" alpha="255" label="0 not assessed"/>
        <paletteEntry value="1" color="#d7191c" alpha="255" label="1 observed event-pair loss"/>
        <paletteEntry value="2" color="#fdae61" alpha="255" label="2 modeled loss confirmed by early POST"/>
        <paletteEntry value="3" color="#7b3294" alpha="255" label="3 post-event decline without event-pair loss"/>
        <paletteEntry value="4" color="#2c7bb6" alpha="255" label="4 no loss evidence"/>
        <paletteEntry value="255" color="#000000" alpha="0" label="255 NoData"/>
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" invertColors="0" colorizeOn="0" colorizeRed="255" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
