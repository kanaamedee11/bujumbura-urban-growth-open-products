// Bujumbura urban growth study (2016 to 2025): Earth Engine export script.
// Paste into https://code.earthengine.google.com , press Run, then press "Run" on the 4 tasks in the Tasks tab.
// No API keys are needed: you are signed in with your own Google account. Output goes to a Drive folder "BujumburaUrban".

var AOI = ee.Geometry.Rectangle([29.28, -3.52, 29.50, -3.28]);   // Bujumbura city and its surroundings (lake included, masked later)
var START_YEAR = 2016, END_YEAR = 2025;
var SCALE = 20;               // metres for the 10-20 m products
var FOLDER = 'BujumburaUrban';
var CRS = 'EPSG:32735';       // UTM 35S

function stackBands(list) {
  return ee.Image(ee.List(list).iterate(function (cur, acc) { return ee.Image(acc).addBands(ee.Image(cur)); }, ee.Image([])));
}
var YEARS = ee.List.sequence(START_YEAR, END_YEAR);
function drySeason(y) {   // June to September, the dry season in Burundi, to limit cloud
  y = ee.Number(y); return [ee.Date.fromYMD(y, 6, 1), ee.Date.fromYMD(y, 10, 1)];
}

// ---------- 1. Dynamic World: dry-season median "built" probability per year, plus the number of observations ----------
var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filterBounds(AOI).select('built');
var built = stackBands(YEARS.map(function (y) {
  var r = drySeason(y);
  return dw.filterDate(r[0], r[1]).median().multiply(100).toUint8().rename(ee.String('built_y').cat(ee.Number(y).format('%d')));
}));
var nobs = stackBands(YEARS.map(function (y) {
  var r = drySeason(y);
  return dw.filterDate(r[0], r[1]).count().min(255).toUint8().rename(ee.String('n_y').cat(ee.Number(y).format('%d')));
}));
Export.image.toDrive({image: built.addBands(nobs).clip(AOI), description: 'bujumbura_dw_built_and_counts',
  folder: FOLDER, region: AOI, scale: SCALE, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'});

// ---------- 2. Land-cover products and terrain ----------
var wc20 = ee.Image('ESA/WorldCover/v100/2020').select('Map').rename('wc2020_class');
var wc21 = ee.Image('ESA/WorldCover/v200/2021').select('Map').rename('wc2021_class');
var srtm = ee.Image('USGS/SRTMGL1_003');
var terrain = ee.Image.cat([wc20, wc21, srtm.rename('elevation_m'), ee.Terrain.slope(srtm).rename('slope_deg')]).toInt16().clip(AOI);
Export.image.toDrive({image: terrain, description: 'bujumbura_landcover_and_terrain',
  folder: FOLDER, region: AOI, scale: SCALE, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'});

// ---------- 3. GHSL built-up surface (100 m, epochs 2000, 2015, 2020) and WorldPop population (100 m, 2015 and 2020) ----------
var ghsl = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_S');
function ghslEpoch(y) { return ghsl.filterDate(y + '-01-01', (y + 1) + '-01-01').first().select('built_surface').rename('ghsl_built_m2_' + y); }
var pop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.eq('country', 'BDI'));
function popYear(y) { return pop.filter(ee.Filter.eq('year', y)).first().select('population').rename('worldpop_' + y); }
var gp = ee.Image.cat([ghslEpoch(2000), ghslEpoch(2015), ghslEpoch(2020), popYear(2015), popYear(2020)]).toFloat().clip(AOI);
Export.image.toDrive({image: gp, description: 'bujumbura_ghsl_worldpop_100m',
  folder: FOLDER, region: AOI, scale: 100, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'});

// ---------- 4. Sentinel-2 true-colour composites for the figures: dry season 2017 and 2025 ----------
function rgb(y) {
  var col = ee.ImageCollection('COPERNICUS/S2_HARMONIZED').filterBounds(AOI)
    .filterDate(y + '-06-01', y + '-10-01').filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .map(function (im) { var qa = im.select('QA60'); return im.updateMask(qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0))); });
  return col.median().select(['B4', 'B3', 'B2']).divide(3000).multiply(255).clamp(0, 255).toUint8()
    .rename(['r' + y, 'g' + y, 'b' + y]);
}
Export.image.toDrive({image: rgb(2017).addBands(rgb(2025)).clip(AOI), description: 'bujumbura_s2_truecolour_2017_2025',
  folder: FOLDER, region: AOI, scale: SCALE, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'});

Map.centerObject(AOI, 11);
Map.addLayer(built.select('built_y2025').selfMask(), {min: 0, max: 100, palette: ['white', 'red']}, 'DW built 2025');
print('Open the Tasks tab (top right) and press Run on the 4 export tasks.');
