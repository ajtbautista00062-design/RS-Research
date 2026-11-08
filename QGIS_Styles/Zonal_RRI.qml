<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="singlebandpseudocolor" band="1" opacity="1" alphaBand="-1" classificationMin="-0.5" classificationMax="1.5" nodataColor="">
      <rasterTransparency/>
      <minMaxOrigin><limits>None</limits><extent>WholeRaster</extent><statAccuracy>Estimated</statAccuracy><cumulativeCutLower>0.02</cumulativeCutLower><cumulativeCutUpper>0.98</cumulativeCutUpper><stdDevFactor>2</stdDevFactor></minMaxOrigin>
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" classificationMode="1" clip="0" minimumValue="-0.5" maximumValue="1.5">
          <item value="-0.5" color="#a50026" alpha="255" label="-0.5 further degradation"/>
          <item value="0" color="#fdae61" alpha="255" label="0 no recovery"/>
          <item value="0.5" color="#ffffbf" alpha="255" label="0.5 half recovered"/>
          <item value="1" color="#66bd63" alpha="255" label="1 back to baseline"/>
          <item value="1.5" color="#006837" alpha="255" label="1.5 above baseline"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" invertColors="0" colorizeOn="0" colorizeRed="255" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
