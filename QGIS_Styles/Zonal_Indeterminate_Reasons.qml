<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <colorPalette>
        <paletteEntry value="0" color="#1a9641" alpha="255" label="0 classified"/>
        <paletteEntry value="1" color="#636363" alpha="255" label="1 no valid observations"/>
        <paletteEntry value="2" color="#fdae61" alpha="255" label="2 insufficient PRE/POST/event support"/>
        <paletteEntry value="3" color="#a50026" alpha="255" label="3 unstable or out-of-range model"/>
        <paletteEntry value="4" color="#7b3294" alpha="255" label="4 baseline fails signal observability gate"/>
        <paletteEntry value="12" color="#f4a582" alpha="255" label="12 insufficient valid zone pixels at every zone size"/>
        <paletteEntry value="255" color="#000000" alpha="0" label="255 NoData"/>
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" invertColors="0" colorizeOn="0" colorizeRed="255" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
