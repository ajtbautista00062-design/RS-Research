<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <colorPalette>
        <paletteEntry value="0" color="#1a9641" alpha="255" label="0 classified"/>
        <paletteEntry value="1" color="#636363" alpha="255" label="1 structurally unavailable"/>
        <paletteEntry value="2" color="#fdae61" alpha="255" label="2 insufficient PRE/POST"/>
        <paletteEntry value="3" color="#a50026" alpha="255" label="3 unstable/nonfinite model"/>
        <paletteEntry value="4" color="#7b3294" alpha="255" label="4 baseline fails configured observability rule"/>
        <paletteEntry value="5" color="#d7191c" alpha="255" label="5 conflicting disturbance"/>
        <paletteEntry value="6" color="#f4a582" alpha="255" label="6 loss near threshold or unsupported"/>
        <paletteEntry value="7" color="#fee090" alpha="255" label="7 late observations unavailable"/>
        <paletteEntry value="8" color="#e66101" alpha="255" label="8 model/observation contradiction"/>
        <paletteEntry value="9" color="#92c5de" alpha="255" label="9 ambiguous recovery"/>
        <paletteEntry value="10" color="#4393c3" alpha="255" label="10 spatial evidence conflict"/>
        <paletteEntry value="11" color="#000000" alpha="0" label="11 outside analysis mask or ancillary exclusion"/>
        <paletteEntry value="255" color="#000000" alpha="0" label="255 NoData"/>
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" invertColors="0" colorizeOn="0" colorizeRed="255" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
